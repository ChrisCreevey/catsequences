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
void clear_memory(void);

char **species_names = '\0', ***sequences = '\0';
int *seqlens = '\0', numspecies = 0, numseqs = 0, num_files=0; 

int main(int argc, char *argv[])
    {
    FILE *infile1 = '\0', *infile2 = '\0', *outfile = '\0';
    char filename[100], string[1000000], string2[1000000], c = '\0';
    int i = 0, j = 0, k=0, found = 0, stopreading = FALSE, aln_len =0, error = FALSE, filenum=0, seqpos=0, previous=FALSE, totlen=0;
    
    if(argc < 2)
        {
        printf("\n\ncatsequences is designed to concatenate sequence alignments. \nUsage: catsequences filename\nwhere filename contains a list of files with sequences to be concatenated\n\n\n");
        exit(1);
        }
    
		
	/* open file with list of files to be concatenated */    
    infile1 = fopen(argv[1], "r");
    
    string[0] ='\0';
    string2[0] = '\0';
    
    /* assign the arry to hold the names of the species (max 1000 in this build) */
	species_names = malloc(10000*sizeof(char*));
    if(species_names == '\0') clean_exit(1);
	
	for(i=0; i<10000; i++)
		{
		species_names[i] = malloc(50*sizeof(char));
		if(species_names[i] == '\0') clean_exit(23);
		species_names[i][0] = '\0';
		}

	
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
	if(seqlens == '\0') clean_exit(6);
	j=0;
	
	while(!feof(infile1) && !error)
		{
        fscanf(infile1, "%s\n", &filename); /* Get the name of the sequene file tobe opened */
        infile2 = fopen(filename, "r"); /* open the sequence file */
        /*printf("%s\n", filename);*/
		
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
					if(c != '\r' && c != '\n' && c != '.' && !feof(infile2))
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
					for(i=0; i<100; i++)   /* Reduce the length of the names to 100 characters */
						{
						if(string[i] == '.')
							{
							i = 100;
							}
						else
							{
							if(string[i] != ' ')
								string2[i] = string[i];
							else
								break;
							}
						}
					string2[i] = '\0';
					strcpy(string, string2);
					
					found = -1;
					for(i=0; i<numspecies; i++)
						{
						if(strcmp(string, species_names[i])==0)
							found = i;
						}
					if(found == -1)
						{
				/*		printf("found species %d: %s\n", numspecies, string); */
						strcpy(species_names[numspecies], string);
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
        fclose(infile2);
		j++;
		}
	fclose(infile1);


	printf("num_species = %d, num_files = %d\n", numspecies, num_files);
	printf("expected concatenated sequence length %d\n", totlen-num_files);


	/*read in the sequences for real */

	/* Assign the array to hold all the sequences */
	sequences = malloc(num_files*sizeof(char**));
    if(sequences == '\0') clean_exit(2);

	for(j=0; j<num_files; j++)
		{
		sequences[j] = malloc(numspecies*sizeof(char*));
		if(sequences[j] == '\0') clean_exit(3);
		for(k=0; k<numspecies; k++) 
			{
			sequences[j][k] = malloc((seqlens[j]+1)*sizeof(char));
			if(sequences[j][k] == '\0') clean_exit(4);
			sequences[j][k][0] = '\0';
			}
		}
		
	printf("Done assignment\n");

	/**** READ IN THE SEQUENCES  ***/


	printf("\n\n\n");
    infile1 = fopen(argv[1], "r");	/* open up the list of files again */
	filenum=0;
	while(!feof(infile1) && !error)
        {		
		fscanf(infile1, "%s\n", &filename); /* read in the name of the file */
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
					string[i] = c;
					if(string[i] != ' ') i++;
					}
				if(i == 1000000) printf("off the end!\n");
				}
			
			string[i] = '\0';
		/*	printf("%s\n", string);*/
			
			
			if(string[0] == '>')
				{
				for(i=0; i<100; i++)   /* Reduce the length of the names to 100 characters */
					{
					if(string[i] == '.')
						{
						string2[i] = '\0';
						i = 100;
						}
					else
						{
						if(string[i] != ' ')
							string2[i] = string[i];
						else
							break;
						}
					}
				string2[i] = '\0';
				strcpy(string, string2);
				
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
					clean_exit(1);
					}
				}
			else
				{
				if(!previous)
					{
					strcat(sequences[filenum][found], string);
					if(strlen(sequences[filenum][found]) > seqlens[filenum]) printf("%d/%d\n", strlen(sequences[filenum][found]), seqlens[filenum]);
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
				if(strlen(sequences[j][i])+1 != seqlens[j]) printf("File %d seqlen = %d. fulllen = %d\n", i, strlen(sequences[j][i]), seqlens[j]);
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
	printf("Finished!\n");
    }
    
void clean_exit(int x)
    {
    printf("Error: out of memory at %d\n", x);
	clear_memory();
    exit(0);
    }
	
void clear_memory(void)
	{
	int i, j;
	if(species_names != '\0')
		{
		for(i=0; i<numspecies; i++)
			{
			if(species_names[i] != '\0')
				{
				free(species_names[i]);
				species_names[i] = '\0';
				}
			}
		free(species_names);
		species_names = '\0';
		}
		
	if(sequences != '\0')
		{
		for(i=0; i<num_files; i++)
			{
			if(sequences[i] != '\0')
				{
				for(j=0; j<numspecies; j++)
					{
					if(sequences[i][j] != '\0')
						{
						free(sequences[i][j]);
						sequences[i][j] = '\0';
						}
					}
				free(sequences[i]);
				sequences[i] = '\0';
				}
			}
		free(sequences);
		sequences = '\0';
		}
	
	if(seqlens != '\0')
		{
		free(seqlens);
		seqlens = '\0';
		}
	}