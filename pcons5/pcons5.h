#ifndef pcons5_h
#define pcons5_h

#define	MAXMETHODS	100		/* Maximum allowable atoms */
#define	MAXFILE		1000	        /* Maximum allowable models */

typedef struct {
  struct
  {
    double cutoffs90_1_5;
    double cutoffs90_3;
    double cutoffs50_1_5;
    double cutoffs50_3;
    char method[200];
  } methods[MAXMETHODS];

} cutoffs;

void readcutoffs(cutoffs *cut,char *file);
//double     LGscore(char *file1,char *file2);
//double*     LGscore_res(char *file1,char *file2);
//double	   LG_pvalue(int N,double MS);
//double	   Levitt_Gerstein(molecule *m1,molecule *m2);
//double	   superimpose_molecules();
//double     LG_pvalueF(int N,double score);//added by Fang
//int        read_molecules_ca();	/* Reads in molecules to be superimposed */
//int        read_molecule(molecule *m);	/* Reads in molecules to be superimposed */
//int        center_molecule(molecule *m);	
//int        copy_matrix(); /* copy matrix f into matrix t */
//int        transpose_matrix(); /* Transpose a 3x3 matrix */
//int        delete_atom(molecule *m1,int num);	
//int        check_molecules(molecule *m1,molecule *m2);
//int        check_molecules2(molecule *m1,molecule *m2);
//int        atom_exist(char *name,int resnum,molecule *m);	
//void       copymolecule(molecule *m1,molecule *m2);


#endif
