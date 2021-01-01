#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TEMP
#define TEMP 2
#endif


void clean_exit(int x);

char **species_names = NULL, ***sequences = NULL;
int *seqlens = NULL, numspecies = 0, numseqs = 0, num_files=0; 

int main(int argc, char *argv[])
    {
    FILE *infile1 = NULL, *infile2 = NULL, *outfile2 = NULL;
    char *filename = NULL, string[1000000], string2[1000000], c = '\0';
    int i = 0, j = 0, k=0, found = 0, aln_len =0, error = FALSE, filenum=0, previous=FALSE, totlen=0;
    
    if(argc < 2)
        {
        printf("\ncatsequences is designed to concatenate sequence alignments. \n\nUsage: catsequences filename\n\tWhere filename contains a list of files with sequences to be concatenated\n\tThe concatenated alignments are outputted to a file called \"allseqs.fas\"\n\tPartition information will be written to the file \"allseqs.partitions.txt\"\n\n\n\tThe sequence names are read up to the first \".\" or \"|\" (or the end of the name)\n\n");
        exit(1);
        }
    
		
	/* open file with list of files to be concatenated */    
	if((infile1 = fopen(argv[1], "r")) == NULL)   /* check to see if the file is there */
	    {                          /* Open the fundamental tree file */
	    fprintf(stderr, "Error: Cannot open list file %s\n", argv[1]);
	    exit(1);
	    }
    
    outfile2 = fopen("allseqs.partitions.txt", "w");
    filename=malloc(1000*sizeof(char));
    filename[0]='\0';
    string[0] ='\0';
    string2[0] = '\0';
    

    /* assign the arry to hold the names of the species (max 10000 in this build) */
	species_names = malloc(10000*sizeof(char*));
    if(species_names == NULL) clean_exit(1);
	
	/* Count the number of file to be read in */
	i=0;
    while(!feof(infile1))
		{
		c= getc(infile1);
		if(c == '\n' || c == '\r') i++;
		}
	rewind(infile1);
	num_files = i;


    /* assign the array to hold the length of the aligned sequences in each file */
	seqlens = malloc(num_files*sizeof(int));
	if(seqlens == NULL) clean_exit(6);
	j=0;
	
	while(!feof(infile1) && !error)
		{
        fscanf(infile1, "%s\n", filename); /* Get the name of the sequene file tobe opened */
		if((infile2 = fopen(filename, "r")) == NULL)   /* check to see if the file is there */
		    {                          /* Open the fundamental tree file */
		    fprintf(stderr, "Error: Cannot open alignment file %s\n", filename);
		    exit(1);
		    }
	   
        /*infile2 = fopen(filename, "r"); */ /* open the sequence file */
        /*printf("%s\n", filename);*/
		
		fprintf(outfile2, "%s\t=\t%d-", filename, totlen-j+1);

		aln_len = -1;
		/* read in the names and see how many species we need for this gene */
		c = getc(infile2);
		while(!feof(infile2) && !error)
            {
			if(c == '>') /* read in the names are */
				{
				i=0;
				string[i] = c;
				i++;
				 while(!feof(infile2) && (c=getc(infile2)) != '\n' && c != '\r' && c != '.')
					{
					if(c != '\r' && c != '\n' && !feof(infile2))
						{
						string[i] = c;
						i++;
						}
					if(i == 1000000) printf("off the end!\n");
					}
				if(c == '.')
					{
					while(!feof(infile2) && (c=getc(infile2)) != '\n' && c != '\r');
					}
				string[i] = '\0';
				
				
				
				if(string[0] == '>')
					{
					string[strcspn(string, " |.")] = '\0'; /* species name is initial string up to [ |.] */
					found = -1;
					for(i=0; i<numspecies; i++)
						{
						if(strcmp(string, species_names[i])==0)
							found = i;
						}
					if(found == -1)
						{
				/*		printf("found species %d: %s\n", numspecies, string); */
						if (!(species_names[numspecies] = strdup(string)))
							{
								perror("strdup");
								exit(EXIT_FAILURE);
							}
						found = numspecies;
						numspecies++;
						}
					}
				}
			i=0;
			
			
			while(!feof(infile2) && (c=getc(infile2)) != '>') /* count the length of the sequences */
				{
				if(c != ' ' && c != '\n' && c != '\r') i++;
				}
			if(!feof(infile2)) i++;
			if(aln_len == -1)
				{
				aln_len = i;
				seqlens[j] = aln_len;
				totlen+=aln_len;
			/*	printf("aln_len = %d seqlens[%d] = %d\n",aln_len, j, seqlens[j]); */
				}
			else
				{
				if(aln_len != i)
					{
					printf("ERROR not all sequences the same length (%d != %d)\n", i , aln_len);
					error=TRUE;
					}
				}
			}
		fprintf(outfile2, "%d;\n", totlen-j-1);
        fclose(infile2);
		j++;
		}
	fclose(infile1);
	fclose(outfile2);

	printf("num_species = %d, num_files = %d\n", numspecies, num_files);
	printf("expected concatenated sequence length %d\n", totlen-num_files);
	printf("\nAll partitions written to \"allseqs.partitions.txt\"\n");


	/*read in the sequences for real */

	/* Assign the array to hold all the sequences */
	sequences = malloc(num_files*sizeof(char**));
    if(sequences == NULL) clean_exit(2);

	for(j=0; j<num_files; j++)
		{
		sequences[j] = malloc(numspecies*sizeof(char*));
		if(sequences[j] == NULL) clean_exit(3);
		for(k=0; k<numspecies; k++) 
			{
			sequences[j][k] = malloc((seqlens[j]+1)*sizeof(char));
			if(sequences[j][k] == NULL) clean_exit(4);
			sequences[j][k][0] = '\0';
			}
		}
		
	printf("Done assignment\n");

	/**** READ IN THE SEQUENCES  ***/


	printf("\n");
    infile1 = fopen(argv[1], "r");	/* open up the list of files again */
	filenum=0;
	while(!feof(infile1) && !error)
        {		
		fscanf(infile1, "%s\n", filename); /* read in the name of the file */
      /*  printf("File opened: %s\n", filename);*/
        infile2 = fopen(filename, "r"); /* open the file */
		previous=FALSE;

		while(!feof(infile2) && !error)
            {
			i=0;
			while(!feof(infile2) && (c=getc(infile2)) != '\n' && c != '\r')
				{
				if(c != '\r' && c != '\n' && !feof(infile2))
					{
					string[i++] = c;
					}
				if(i == 1000000) printf("off the end!\n");
				}
			
			string[i] = '\0';
		/*	printf("%s\n", string);*/
			
			
			if(string[0] == '>')
				{
				string[strcspn(string, " |.")] = '\0'; /* species name is initial string up to [ |.] */
				found = -1;
				for(i=0; i<numspecies; i++)
					{
					if(strcmp(string, species_names[i]) == 0)
						{
						found = i;
						if(strlen(sequences[filenum][found]) > 0)
							{
							previous=TRUE;
							printf("Found an extra copy of %s in file %s.... ignoring this copy\n", species_names[i], filename);
							}
						i=numspecies;
						}
					}
				if(found == -1)
					{
					printf("error: %s not found\n", string);
					exit(EXIT_FAILURE);
					}
				}
			else
				{
				if(!previous)
					{
					strcat(sequences[filenum][found], string);
					if(strlen(sequences[filenum][found]) > seqlens[filenum]) printf("%lu/%d\n", strlen(sequences[filenum][found]), seqlens[filenum]);
					}
				}


			i=0;
			}
        fclose(infile2);
		filenum++;
        }
    fclose(infile1);



/*	for(i=0; i<filenum; i++)
		{
		printf("File %d seqlen = %d. fulllen = %d\n", i, strlen(sequences[i][0]), seqlens[i]);
		}
*/
    infile1 = fopen("allseqs.fas", "w");
    
    found = 0;
    for(i=0; i<numspecies; i++)
        {
        fprintf(infile1, "%s\n", species_names[i]);
        for(j=0; j<num_files; j++)
            {
            if(strlen(sequences[j][i]) != 0)
                {
                fprintf(infile1, "%s", sequences[j][i]);
				if(strlen(sequences[j][i])+1 != seqlens[j]) printf("File %d seqlen = %lu. fulllen = %d\n", i, strlen(sequences[j][i]), seqlens[j]);
/*				printf("%d/%d\t", strlen(sequences[j][i]), seqlens[j]);*/
                }
            else
                {
/*				printf("-%d\t", seqlens[j]-1); */
                for(k=0; k<(seqlens[j]-1); k++)
                    {
                    fprintf(infile1, "?");
                    }
                }
            
            
            }
        fprintf(infile1, "\n");
	/*	printf("\n"); */
        }
    fclose(infile1);
    free(filename);
	printf("Finished! Cancatenated alignment written to \"allseqs.fas\"\n\n");
    }
    
void clean_exit(int x)
    {
    printf("Error: out of memory at %d\n", x);
    exit(EXIT_FAILURE);
    }
