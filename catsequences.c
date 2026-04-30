/*
 * catsequences - concatenate multiple FASTA alignments for supermatrix
 * phylogenetic analyses.
 *
 * Usage:
 *   catsequences [OPTIONS] LISTFILE
 *   catsequences [OPTIONS] FILE1.fasta FILE2.fasta ...
 *
 * The list file contains one alignment filename per line (glob wildcards
 * are supported).  Alternatively pass the FASTA files directly on the
 * command line.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <glob.h>
#include <errno.h>

#ifndef TRUE
#  define TRUE  1
#  define FALSE 0
#endif

/* Growing-array initial capacities */
#define INIT_SPECIES 256
#define INIT_FILES   256
#define INIT_BUFLEN  4096

/* Partition output formats */
#define PFMT_GENERIC 0
#define PFMT_RAXML   1
#define PFMT_IQTREE  2
#define PFMT_NEXUS   3

/* Alignment output formats */
#define AFMT_FASTA   0
#define AFMT_NEXUS   1
#define AFMT_PHYLIP  2

/* Duplicate-handling modes */
#define DUP_WARN     0
#define DUP_ERROR    1
#define DUP_KEEPLAST 2

/* ── global arrays (freed at exit) ───────────────────────────── */
static char  **species_names      = NULL;
static int     numspecies         = 0;
static int     species_cap        = 0;

static char  **filelist           = NULL;
static int     num_files          = 0;
static int     files_cap          = 0;

static int    *seqlens            = NULL; /* actual sequence length per file */
static int    *file_active        = NULL; /* 1 if file passes filters        */
static int    *file_species_count = NULL; /* species found per file          */

static char ***sequences          = NULL; /* sequences[file][species]        */

/* ── runtime options ─────────────────────────────────────────── */
static const char *name_delimiters = " |.";
static char        missing_char    = '?';
static const char *out_prefix      = "allseqs";
static int         quiet           = FALSE;
static int         verbose         = FALSE;
static int         part_format     = PFMT_GENERIC;
static int         aln_format      = AFMT_FASTA;
static int         dup_mode        = DUP_WARN;
static double      min_taxa_frac   = 0.0;
static const char *include_file    = NULL;
static const char *exclude_file    = NULL;
static int         write_occ_csv   = FALSE;

static char **include_list = NULL;
static int    include_count = 0;
static char **exclude_list = NULL;
static int    exclude_count = 0;

/* ═══════════════════════ helpers ════════════════════════════════ */

static void die(const char *msg)
{
    fprintf(stderr, "Error: %s\n", msg);
    exit(EXIT_FAILURE);
}

/*
 * Read one line from fp into *buf (grows via realloc as needed).
 * Strips the trailing \r and/or \n.
 * Returns the number of characters stored (>=0), or -1 on EOF with no data.
 */
static int read_line(FILE *fp, char **buf, size_t *cap)
{
    size_t len = 0;
    int c;

    if (*buf == NULL) {
        *cap = INIT_BUFLEN;
        *buf = malloc(*cap);
        if (!*buf) die("out of memory");
    }

    for (;;) {
        c = getc(fp);
        if (c == EOF) {
            if (len == 0) return -1;
            break;
        }
        if (c == '\n') break;
        if (c == '\r') {
            int nxt = getc(fp);
            if (nxt != '\n' && nxt != EOF) ungetc(nxt, fp);
            break;
        }
        if (len + 1 >= *cap) {
            *cap *= 2;
            char *tmp = realloc(*buf, *cap);
            if (!tmp) die("out of memory");
            *buf = tmp;
        }
        (*buf)[len++] = (char)c;
    }
    (*buf)[len] = '\0';
    return (int)len;
}

/* Return basename of path with the extension stripped. */
static void gene_name(const char *path, char *out, size_t outsz)
{
    const char *base = strrchr(path, '/');
    base = base ? base + 1 : path;
    const char *dot = strrchr(base, '.');
    size_t len = dot ? (size_t)(dot - base) : strlen(base);
    if (len >= outsz) len = outsz - 1;
    memcpy(out, base, len);
    out[len] = '\0';
}

/* Load one name per line from a file into *list / *count. */
static void load_name_list(const char *path, char ***list, int *count)
{
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open %s: %s\n", path, strerror(errno));
        exit(EXIT_FAILURE);
    }
    int cap = 64;
    *list = malloc((size_t)cap * sizeof(char *));
    if (!*list) die("out of memory");
    *count = 0;

    char *buf = NULL;
    size_t bufcap = 0;
    while (read_line(fp, &buf, &bufcap) >= 0) {
        if (buf[0] == '\0') continue;
        if (*count >= cap) {
            cap *= 2;
            char **tmp = realloc(*list, (size_t)cap * sizeof(char *));
            if (!tmp) die("out of memory");
            *list = tmp;
        }
        (*list)[*count] = strdup(buf);
        if (!(*list)[*count]) die("out of memory");
        (*count)++;
    }
    free(buf);
    fclose(fp);
}

/*
 * Return TRUE if name (with or without leading '>') matches a list entry.
 * List entries are stored without '>'.
 */
static int in_list(const char *name, char **list, int count)
{
    const char *nm = (name[0] == '>') ? name + 1 : name;
    for (int i = 0; i < count; i++)
        if (strcmp(nm, list[i]) == 0) return TRUE;
    return FALSE;
}

/* Append one path to the global filelist, growing as needed. */
static void add_file(const char *path)
{
    if (num_files >= files_cap) {
        files_cap = files_cap ? files_cap * 2 : INIT_FILES;
        char **tmp = realloc(filelist, (size_t)files_cap * sizeof(char *));
        if (!tmp) die("out of memory");
        filelist = tmp;
    }
    filelist[num_files] = strdup(path);
    if (!filelist[num_files]) die("out of memory");
    num_files++;
}

/* Expand a glob pattern and add all matches; treat as literal if no match. */
static void add_pattern(const char *pattern)
{
    glob_t g;
    int r = glob(pattern, GLOB_TILDE, NULL, &g);
    if (r == 0) {
        for (size_t i = 0; i < g.gl_pathc; i++)
            add_file(g.gl_pathv[i]);
        globfree(&g);
    } else {
        if (r != GLOB_NOMATCH) globfree(&g);
        add_file(pattern);
    }
}

/* Read a list file and expand each line as a glob pattern. */
static void load_file_list(const char *listpath)
{
    FILE *fp = fopen(listpath, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open list file %s: %s\n",
                listpath, strerror(errno));
        exit(EXIT_FAILURE);
    }
    char *buf = NULL;
    size_t cap = 0;
    while (read_line(fp, &buf, &cap) >= 0) {
        if (buf[0] == '\0') continue;
        add_pattern(buf);
    }
    free(buf);
    fclose(fp);
}

/* Pre-flight: verify every file is readable; report all missing at once. */
static void preflight_check(void)
{
    int bad = FALSE;
    for (int i = 0; i < num_files; i++) {
        FILE *fp = fopen(filelist[i], "r");
        if (!fp) {
            fprintf(stderr, "Error: Cannot open %s: %s\n",
                    filelist[i], strerror(errno));
            bad = TRUE;
        } else {
            fclose(fp);
        }
    }
    if (bad) exit(EXIT_FAILURE);
}

/* Find species by name, or add a new entry.  Returns the index. */
static int find_or_add_species(const char *name)
{
    for (int i = 0; i < numspecies; i++)
        if (strcmp(name, species_names[i]) == 0) return i;

    if (numspecies >= species_cap) {
        species_cap = species_cap ? species_cap * 2 : INIT_SPECIES;
        char **tmp = realloc(species_names,
                             (size_t)species_cap * sizeof(char *));
        if (!tmp) die("out of memory");
        species_names = tmp;
    }
    species_names[numspecies] = strdup(name);
    if (!species_names[numspecies]) die("out of memory");
    return numspecies++;
}

/* ═════════════════════════ first pass ═══════════════════════════
 * Scan each FASTA file to build the species list and measure
 * sequence lengths.  No sequences are stored yet.
 * ═══════════════════════════════════════════════════════════════ */
static void first_pass(void)
{
    seqlens            = calloc((size_t)num_files, sizeof(int));
    file_active        = malloc((size_t)num_files * sizeof(int));
    file_species_count = calloc((size_t)num_files, sizeof(int));
    if (!seqlens || !file_active || !file_species_count) die("out of memory");
    for (int i = 0; i < num_files; i++) file_active[i] = TRUE;

    char *buf = NULL;
    size_t cap = 0;

    for (int j = 0; j < num_files; j++) {
        FILE *fp = fopen(filelist[j], "r");
        if (!fp) {
            fprintf(stderr, "Error: Cannot open %s\n", filelist[j]);
            exit(EXIT_FAILURE);
        }

        int aln_len = -1;
        int in_seq  = FALSE;
        int seq_len = 0;

        while (read_line(fp, &buf, &cap) >= 0) {
            if (buf[0] == '>') {
                /* Commit length of the sequence just finished */
                if (in_seq) {
                    if (aln_len == -1) {
                        aln_len    = seq_len;
                        seqlens[j] = aln_len;
                    } else if (seq_len != aln_len) {
                        fprintf(stderr,
                            "ERROR: sequences differ in length in %s"
                            " (%d vs %d)\n",
                            filelist[j], seq_len, aln_len);
                        exit(EXIT_FAILURE);
                    }
                    seq_len = 0;
                }
                in_seq = TRUE;

                /* Truncate the header at the first delimiter (after '>') */
                buf[1 + strcspn(buf + 1, name_delimiters)] = '\0';
                find_or_add_species(buf);
                file_species_count[j]++;

            } else if (in_seq) {
                /* Count non-whitespace characters in sequence lines */
                for (const char *p = buf; *p; p++)
                    if (*p != ' ' && *p != '\t') seq_len++;
            }
        }

        /* Commit the last sequence in the file */
        if (in_seq) {
            if (aln_len == -1) {
                aln_len    = seq_len;
                seqlens[j] = aln_len;
            } else if (seq_len != aln_len) {
                fprintf(stderr,
                    "ERROR: sequences differ in length in %s (%d vs %d)\n",
                    filelist[j], seq_len, aln_len);
                exit(EXIT_FAILURE);
            }
        }

        fclose(fp);

        if (verbose)
            printf("Scanned %s: %d bp, %d species\n",
                   filelist[j], seqlens[j], file_species_count[j]);
    }
    free(buf);

    /* Apply --min-taxa filter now that we know the total species count */
    if (min_taxa_frac > 0.0 && numspecies > 0) {
        for (int j = 0; j < num_files; j++) {
            double frac = (double)file_species_count[j] / (double)numspecies;
            if (frac < min_taxa_frac) {
                file_active[j] = FALSE;
                if (!quiet)
                    printf("Skipping %s: %.1f%% species occupancy"
                           " < required %.1f%%\n",
                           filelist[j], frac * 100.0,
                           min_taxa_frac * 100.0);
            }
        }
    }
}

/* ════════════════════════ second pass ═══════════════════════════
 * Read the actual sequence data into memory.
 * ═══════════════════════════════════════════════════════════════ */
static void second_pass(void)
{
    sequences = calloc((size_t)num_files, sizeof(char **));
    if (!sequences) die("out of memory");

    for (int j = 0; j < num_files; j++) {
        sequences[j] = calloc((size_t)numspecies, sizeof(char *));
        if (!sequences[j]) die("out of memory");
        for (int k = 0; k < numspecies; k++) {
            sequences[j][k] = malloc((size_t)(seqlens[j] + 1));
            if (!sequences[j][k]) die("out of memory");
            sequences[j][k][0] = '\0';
        }
    }

    char *buf = NULL;
    size_t cap = 0;

    for (int j = 0; j < num_files; j++) {
        FILE *fp = fopen(filelist[j], "r");
        if (!fp) {
            fprintf(stderr, "Error: Cannot reopen %s\n", filelist[j]);
            exit(EXIT_FAILURE);
        }

        int found    = -1;
        int previous = FALSE;

        while (read_line(fp, &buf, &cap) >= 0) {
            if (buf[0] == '>') {
                buf[1 + strcspn(buf + 1, name_delimiters)] = '\0';

                found    = -1;
                previous = FALSE;
                for (int i = 0; i < numspecies; i++) {
                    if (strcmp(buf, species_names[i]) == 0) {
                        found = i;
                        break;
                    }
                }
                if (found < 0) {
                    fprintf(stderr,
                        "Error: species '%s' not in master list\n", buf);
                    exit(EXIT_FAILURE);
                }

                /* Handle duplicates */
                if (sequences[j][found][0] != '\0') {
                    switch (dup_mode) {
                    case DUP_WARN:
                        printf("Duplicate %s in %s: ignoring extra copy\n",
                               species_names[found] + 1, filelist[j]);
                        previous = TRUE;
                        break;
                    case DUP_ERROR:
                        fprintf(stderr,
                            "Error: duplicate species %s in %s\n",
                            species_names[found] + 1, filelist[j]);
                        exit(EXIT_FAILURE);
                    case DUP_KEEPLAST:
                        sequences[j][found][0] = '\0';
                        break;
                    }
                }

            } else if (found >= 0 && !previous) {
                /* Append sequence line, guarding against overflow */
                size_t have = strlen(sequences[j][found]);
                size_t room = (size_t)seqlens[j] - have;
                strncat(sequences[j][found], buf, room);
            }
        }
        fclose(fp);

        if (verbose)
            printf("Read %s\n", filelist[j]);
    }
    free(buf);
}

/* ════════════════════════ output helpers ════════════════════════ */

static int count_active(void)
{
    int n = 0;
    for (int j = 0; j < num_files; j++) if (file_active[j]) n++;
    return n;
}

static int total_aln_len(void)
{
    int len = 0;
    for (int j = 0; j < num_files; j++)
        if (file_active[j]) len += seqlens[j];
    return len;
}

/* Return TRUE if species i should appear in the output. */
static int species_shown(int i)
{
    if (include_count > 0 &&
        !in_list(species_names[i], include_list, include_count))
        return FALSE;
    if (exclude_count > 0 &&
         in_list(species_names[i], exclude_list, exclude_count))
        return FALSE;
    return TRUE;
}

static int count_shown_species(void)
{
    int n = 0;
    for (int i = 0; i < numspecies; i++) if (species_shown(i)) n++;
    return n;
}

/* Longest displayed name (without '>') among shown species. */
static int max_displayed_name_len(void)
{
    int mx = 0;
    for (int i = 0; i < numspecies; i++) {
        if (!species_shown(i)) continue;
        int l = (int)strlen(species_names[i]) - 1; /* skip '>' */
        if (l > mx) mx = l;
    }
    return mx;
}

/* ════════════════════════ write partitions ══════════════════════ */
static void write_partitions(void)
{
    char path[4096];
    snprintf(path, sizeof(path), "%s.partitions.txt", out_prefix);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot create %s\n", path);
        exit(EXIT_FAILURE);
    }

    int  totlen = 0;
    char gname[512];

    if (part_format == PFMT_NEXUS)
        fprintf(fp, "#NEXUS\nbegin sets;\n");

    for (int j = 0; j < num_files; j++) {
        if (!file_active[j]) continue;
        int start = totlen + 1;
        totlen   += seqlens[j];
        gene_name(filelist[j], gname, sizeof(gname));

        switch (part_format) {
        case PFMT_GENERIC:
            fprintf(fp, "%s\t=\t%d-%d;\n", filelist[j], start, totlen);
            break;
        case PFMT_RAXML:
            fprintf(fp, "DNA, %s = %d-%d\n", gname, start, totlen);
            break;
        case PFMT_IQTREE:
            fprintf(fp, "%s = %d-%d\n", gname, start, totlen);
            break;
        case PFMT_NEXUS:
            fprintf(fp, "\tcharset %s = %d-%d;\n", gname, start, totlen);
            break;
        }
    }

    if (part_format == PFMT_NEXUS) fprintf(fp, "end;\n");
    fclose(fp);
    if (!quiet) printf("Partitions written to \"%s\"\n", path);
}

/* ════════════════════════ write alignment ═══════════════════════ */
static void write_alignment(void)
{
    const char *ext = aln_format == AFMT_NEXUS  ? "nex" :
                      aln_format == AFMT_PHYLIP  ? "phy" : "fas";
    char path[4096];
    snprintf(path, sizeof(path), "%s.%s", out_prefix, ext);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot create %s\n", path);
        exit(EXIT_FAILURE);
    }

    int nsp     = count_shown_species();
    int tot     = total_aln_len();
    int namelen = max_displayed_name_len();

    /* Format-specific headers */
    if (aln_format == AFMT_PHYLIP) {
        fprintf(fp, " %d %d\n", nsp, tot);
    } else if (aln_format == AFMT_NEXUS) {
        fprintf(fp, "#NEXUS\nbegin data;\n");
        fprintf(fp, "\tdimensions ntax=%d nchar=%d;\n", nsp, tot);
        fprintf(fp, "\tformat datatype=dna missing=%c gap=-;\n", missing_char);
        fprintf(fp, "matrix\n");
    }

    for (int i = 0; i < numspecies; i++) {
        if (!species_shown(i)) continue;

        const char *dname = species_names[i] + 1; /* skip '>' */

        if (aln_format == AFMT_FASTA) {
            fprintf(fp, "%s\n", species_names[i]); /* include '>' */
        } else if (aln_format == AFMT_PHYLIP) {
            fprintf(fp, "%-*s ", namelen, dname);
        } else { /* NEXUS */
            fprintf(fp, "%-*s\t", namelen, dname);
        }

        for (int j = 0; j < num_files; j++) {
            if (!file_active[j]) continue;
            if (sequences[j][i][0] != '\0') {
                fprintf(fp, "%s", sequences[j][i]);
            } else {
                for (int k = 0; k < seqlens[j]; k++)
                    fputc(missing_char, fp);
            }
        }
        fputc('\n', fp);
    }

    if (aln_format == AFMT_NEXUS) fprintf(fp, ";\nend;\n");
    fclose(fp);
    if (!quiet) printf("Alignment written to \"%s\"\n", path);
}

/* ════════════════════════ write occupancy ═══════════════════════ */
static void write_occupancy(void)
{
    char path[4096];
    snprintf(path, sizeof(path), "%s.occupancy.txt", out_prefix);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "Warning: cannot create %s\n", path);
        return;
    }

    int  active       = count_active();
    int  tot          = total_aln_len();
    long total_cells  = 0;
    long present_cells = 0;

    fprintf(fp, "Species\tPresent_genes\tMissing_chars\tMissing_pct\n");

    for (int i = 0; i < numspecies; i++) {
        if (!species_shown(i)) continue;
        int present = 0, miss_chars = 0;
        for (int j = 0; j < num_files; j++) {
            if (!file_active[j]) continue;
            total_cells++;
            if (sequences[j][i][0] != '\0') {
                present++;
                present_cells++;
            } else {
                miss_chars += seqlens[j];
            }
        }
        double pct = tot > 0 ? 100.0 * miss_chars / tot : 0.0;
        fprintf(fp, "%s\t%d/%d\t%d\t%.1f%%\n",
                species_names[i] + 1, present, active, miss_chars, pct);
    }

    if (total_cells > 0)
        fprintf(fp,
            "\nMatrix completeness: %.1f%%"
            " (%ld/%ld gene-species cells present)\n",
            100.0 * (double)present_cells / (double)total_cells,
            present_cells, total_cells);

    fclose(fp);
    if (!quiet) printf("Occupancy summary written to \"%s\"\n", path);
}

static void write_occupancy_csv(void)
{
    char path[4096];
    snprintf(path, sizeof(path), "%s.occupancy.csv", out_prefix);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "Warning: cannot create %s\n", path);
        return;
    }

    char gname[512];
    fprintf(fp, "species");
    for (int j = 0; j < num_files; j++) {
        if (!file_active[j]) continue;
        gene_name(filelist[j], gname, sizeof(gname));
        fprintf(fp, ",%s", gname);
    }
    fputc('\n', fp);

    for (int i = 0; i < numspecies; i++) {
        if (!species_shown(i)) continue;
        fprintf(fp, "%s", species_names[i] + 1);
        for (int j = 0; j < num_files; j++) {
            if (!file_active[j]) continue;
            fprintf(fp, ",%d", sequences[j][i][0] != '\0' ? 1 : 0);
        }
        fputc('\n', fp);
    }

    fclose(fp);
    if (!quiet) printf("Occupancy CSV written to \"%s\"\n", path);
}

/* ════════════════════════ cleanup ════════════════════════════════ */
static void free_all(void)
{
    if (species_names) {
        for (int i = 0; i < numspecies; i++) free(species_names[i]);
        free(species_names);
    }
    if (sequences) {
        for (int j = 0; j < num_files; j++) {
            if (sequences[j]) {
                for (int k = 0; k < numspecies; k++)
                    free(sequences[j][k]);
                free(sequences[j]);
            }
        }
        free(sequences);
    }
    if (filelist) {
        for (int i = 0; i < num_files; i++) free(filelist[i]);
        free(filelist);
    }
    free(seqlens);
    free(file_active);
    free(file_species_count);
    if (include_list) {
        for (int i = 0; i < include_count; i++) free(include_list[i]);
        free(include_list);
    }
    if (exclude_list) {
        for (int i = 0; i < exclude_count; i++) free(exclude_list[i]);
        free(exclude_list);
    }
}

/* ════════════════════════ usage ══════════════════════════════════ */
static void usage(const char *prog)
{
    printf(
        "\ncatsequences - concatenate FASTA alignments for supermatrix"
        " phylogenetics\n"
        "\nUsage:\n"
        "  %s [OPTIONS] LISTFILE\n"
        "  %s [OPTIONS] FILE1.fasta FILE2.fasta ...\n"
        "\nOptions:\n"
        "  -d DELIMS          Name delimiter characters"
        " (default: \" |.\")\n"
        "  -o PREFIX          Output filename prefix (default: allseqs)\n"
        "  -m CHAR            Missing-data character (default: ?)\n"
        "  -q                 Quiet: suppress progress output\n"
        "  -v                 Verbose: print per-file progress\n"
        "  -f FORMAT          Partition format:\n"
        "                       generic (default) | raxml | iqtree | nexus\n"
        "  -F FORMAT          Alignment format:\n"
        "                       fasta (default) | nexus | phylip\n"
        "  -O                 Write occupancy matrix CSV\n"
        "  --min-taxa N       Minimum fraction of species required per"
        " gene [0.0-1.0]\n"
        "  --on-duplicate M   Duplicate-entry handling:\n"
        "                       warn (default) | error | keep-last\n"
        "  --include FILE     Only include species listed in FILE\n"
        "  --exclude FILE     Exclude species listed in FILE\n"
        "\nOutputs:\n"
        "  <PREFIX>.fas / .nex / .phy   Concatenated alignment\n"
        "  <PREFIX>.partitions.txt      Partition information\n"
        "  <PREFIX>.occupancy.txt       Per-species occupancy summary\n"
        "  <PREFIX>.occupancy.csv       Occupancy matrix (with -O)\n"
        "\n",
        prog, prog);
}

/* ═════════════════════════ main ══════════════════════════════════ */
int main(int argc, char *argv[])
{
    static const struct option long_opts[] = {
        {"min-taxa",     required_argument, NULL, 1},
        {"on-duplicate", required_argument, NULL, 2},
        {"include",      required_argument, NULL, 3},
        {"exclude",      required_argument, NULL, 4},
        {"help",         no_argument,       NULL, 'h'},
        {NULL, 0, NULL, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "d:o:m:qvf:F:Oh",
                              long_opts, NULL)) != -1) {
        switch (opt) {
        case 'd': name_delimiters = optarg;    break;
        case 'o': out_prefix      = optarg;    break;
        case 'm': missing_char    = optarg[0]; break;
        case 'q': quiet           = TRUE;      break;
        case 'v': verbose         = TRUE;      break;
        case 'O': write_occ_csv   = TRUE;      break;
        case 'f':
            if      (strcmp(optarg, "generic") == 0) part_format = PFMT_GENERIC;
            else if (strcmp(optarg, "raxml")   == 0) part_format = PFMT_RAXML;
            else if (strcmp(optarg, "iqtree")  == 0) part_format = PFMT_IQTREE;
            else if (strcmp(optarg, "nexus")   == 0) part_format = PFMT_NEXUS;
            else {
                fprintf(stderr, "Unknown partition format: %s\n", optarg);
                exit(EXIT_FAILURE);
            }
            break;
        case 'F':
            if      (strcmp(optarg, "fasta")  == 0) aln_format = AFMT_FASTA;
            else if (strcmp(optarg, "nexus")  == 0) aln_format = AFMT_NEXUS;
            else if (strcmp(optarg, "phylip") == 0) aln_format = AFMT_PHYLIP;
            else {
                fprintf(stderr, "Unknown alignment format: %s\n", optarg);
                exit(EXIT_FAILURE);
            }
            break;
        case 1: /* --min-taxa */
            min_taxa_frac = atof(optarg);
            if (min_taxa_frac < 0.0 || min_taxa_frac > 1.0) {
                fprintf(stderr,
                    "Error: --min-taxa must be in [0.0, 1.0]\n");
                exit(EXIT_FAILURE);
            }
            break;
        case 2: /* --on-duplicate */
            if      (strcmp(optarg, "warn")      == 0) dup_mode = DUP_WARN;
            else if (strcmp(optarg, "error")     == 0) dup_mode = DUP_ERROR;
            else if (strcmp(optarg, "keep-last") == 0) dup_mode = DUP_KEEPLAST;
            else {
                fprintf(stderr, "Unknown --on-duplicate mode: %s\n", optarg);
                exit(EXIT_FAILURE);
            }
            break;
        case 3: include_file = optarg; break;
        case 4: exclude_file = optarg; break;
        case 'h': usage(argv[0]); return EXIT_SUCCESS;
        default:  usage(argv[0]); exit(EXIT_FAILURE);
        }
    }

    if (optind >= argc) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    /* Load include / exclude species lists if provided */
    if (include_file)
        load_name_list(include_file, &include_list, &include_count);
    if (exclude_file)
        load_name_list(exclude_file, &exclude_list, &exclude_count);

    /*
     * Determine whether the remaining argument(s) are a list file or
     * direct FASTA files.  A single argument whose first byte is '>'
     * is treated as a FASTA file; otherwise it is a list file.
     * Multiple arguments are always treated as direct FASTA files.
     */
    if (argc - optind == 1) {
        FILE *probe = fopen(argv[optind], "r");
        if (!probe) {
            fprintf(stderr, "Error: Cannot open %s: %s\n",
                    argv[optind], strerror(errno));
            exit(EXIT_FAILURE);
        }
        int first = getc(probe);
        fclose(probe);
        if (first == '>') {
            add_file(argv[optind]);       /* single FASTA file */
        } else {
            load_file_list(argv[optind]); /* list file         */
        }
    } else {
        /* Multiple positional arguments → direct FASTA files */
        for (int i = optind; i < argc; i++)
            add_pattern(argv[i]);
    }

    if (num_files == 0) {
        fprintf(stderr, "Error: No input files specified\n");
        exit(EXIT_FAILURE);
    }

    preflight_check();

    if (!quiet)
        printf("catsequences: processing %d file(s)\n", num_files);

    first_pass();

    if (!quiet)
        printf("Found %d species across %d file(s)\n", numspecies, num_files);

    second_pass();

    write_partitions();
    write_alignment();
    write_occupancy();
    if (write_occ_csv) write_occupancy_csv();

    if (!quiet) {
        printf("Total concatenated alignment length: %d bp\n",
               total_aln_len());
        printf("Done.\n");
    }

    free_all();
    return EXIT_SUCCESS;
}
