#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "molecule.h"
#include "nets.h"
int get_residue_type623(char* res,char ss,double surface,double buried_cutoff,double exposed_cutoff);
int get_residue_type622(char* res,char ss,double surface,double buried_cutoff);
int get_residue_type632(char* res,char ss,double surface,double buried_cutoff);

main(argc,argv)		/* Main routine */
     int argc;
     char *argv[];
{
  molecule      m[1];
  network       net[5];
  int        windowsize=9;
  int        total_res_contacts_count=0;
  int tmp=0;
  int           surface_total=0;
  int           restotal=0;
  int           res_in_surffile=0;
  int           residues=0;
  int        rescount=0;
  int           res_contact_map[MAXRES][MAXRES]={{}};
  int           res[6][6]={{}};
  int           i,j,k,l,n;
  int           first=0;
  int           current_res_j=0;
  int           current_res_i=0;
  int           current_atomtype_j=0;
  int           current_atomtype_i=0;
  int           current_atomtype3_j=0;
  int           current_atomtype3_i=0;
  int           atomtotal=0;
  int           atoms=0;
  int           atomcount=0;
  int           atomcount3=0;
  int           atomcontacts[MAXRES][14][14]={{{}}};
  int           atomcontacts3[MAXRES][3][3]={{{}}};
  int           atom[14][14]={{}};
  int           atom3[3][3]={{}};
  int           alternative_ss_calc=0;
  int           check=0;
  int           residuetypes=6;
  int           atomtypes=13;
  int           noprof=1;
  int           restype1=0;
  int           restype2=0;
  int           max_j,max_k,max_l;
  FILE          *fp;
  double        total_res_contacts=0;
  double temp;
  double           surface_total_prof=0;
  double        surfaces[MAXRES]={};
  double        surface75_prof[20]={};
  double        surface75[6]={};
  double        surface50_prof[20]={};
  double        surface50[6]={};
  double        surface25_prof[20]={};
  double        surface25[6]={};
  double        surface100_prof[20]={};
  double        surface100[6]={};
  double     ss=0;
  double        sheet[MAXRES];
  double        resfrac[666]={};
  double        resfrac6[21]={};
  //  double        resfrac36[666]={};
  double        rescontacts_prof[36][36]={{}};
  double        rescontacts_prof6[6][6]={{}};
  double        rescontacts_prof622[24][24]={{}};
  double        rescontacts_prof623[36][36]={{}};
  double        rescontacts_prof632[36][36]={{}};
  double        prof[MAXRES][22]={{}};
  double        mean=0;
  double        helix[MAXRES];
  double        dist;
  double        cutoff=25; //residue cutoff 7
  double        coil[MAXRES];
  double     atomfrac[91]={};
  double     atomfrac3[6]={};
  double        atomcutoff=25; // atom cutoff 5;
  double        exposed_cutoff=70;
  double        buried_cutoff=30;
  double        tmp_sum=0;
  double        parameters[137];
  double        predicted_quality=0;
  double        quality[5];
  char          temp_surface[10];
  char          res1[2];
  char          res2[2];
  char          surfacefile[1500]="undef";
  char          stridefile[1500]="undef";
  char *stride;
  char          psipredfile[1500]="undef";
  char *psipred;
  char          profseq[MAXRES];
  char          proffile[1500]="undef";
  char          line_flag[10];
  char          buff[1000];
  char          basename[1000]="undef";
  char          tmp_str[100000];
  char          tmp_short_str[100];
  char          net1[200],net2[200],net3[200],net4[200],net5[200];
  //char          netpath[2000];
  //float        sum=0;
  /* Parse command line for PDB filename */



  i=1;
  //strcpy(net_path,argv[0]);
  while(i<argc)
    {
      //printf("%d %d %s\n",i,argc,argv[i]);
      if (strcmp(argv[i],"-pdb")==0)
	{
	  strcpy(m[0].filename,argv[i+1]);
	  i++;
	}
      else if (strcmp(argv[i],"-surf")==0)
	{  
	  strcpy(surfacefile,argv[i+1]);
	  i++;
	}
      else if (strcmp(argv[i],"-stride")==0)
	{  
	  strcpy(stridefile,argv[i+1]);
	  i++;
	}
      else if (strcmp(argv[i],"-psipred")==0)
	{  
	  strcpy(psipredfile,argv[i+1]);
	  i++;
	}
      else if (strcmp(argv[i],"-noprof")==0)
	{
	  noprof=1;
	}
      else if (strcmp(argv[i],"-prof")==0)
	{
	  strcpy(proffile,argv[i+1]);
	  i++;
	}
      else if (strcmp(argv[i],"-w")==0)
	{  
	  windowsize=atoi(argv[i+1]);
	  i++;
	}
      else if (strcmp(argv[i],"-rc")==0)
	{  
	  cutoff=strtod(argv[i+1],NULL); //argv[i+1]+strlen(argv[i+1]));
	  cutoff=cutoff*cutoff;
	  i++;
	}
      else if (strcmp(argv[i],"-ac")==0)
	{  
	  atomcutoff=strtod(argv[i+1],NULL); //argv[i+1]+strlen(argv[i+1]));
	  atomcutoff=atomcutoff*atomcutoff;
	  i++;
	}
      else if (strcmp(argv[i],"-ec")==0)
	{  
	  exposed_cutoff=strtod(argv[i+1],NULL); //argv[i+1]+strlen(argv[i+1]));
	  i++;
	}
      else if (strcmp(argv[i],"-bc")==0)
	{  
	  buried_cutoff=strtod(argv[i+1],NULL); //argv[i+1]+strlen(argv[i+1]));
	  i++;
	}
      else if (strcmp(argv[i],"-t")==0)
	{  
	  residuetypes=atoi(argv[i+1]);
	  if(residuetypes != 6 && residuetypes != 20 && residuetypes != 623 && residuetypes != 632 && residuetypes != 622)
	    {
	      fprintf(stderr,"%d is an unallowed number of residues types ...\nOnly 6, 20, 623, 622 or 632 residue types are allowed!\n",residuetypes);
	      exit(1);
	    }
	  i++;
	}
       else if (strcmp(argv[i],"-atomtype")==0)
	{  
	  atomtypes=atoi(argv[i+1]);
	  if(atomtypes != 3 && atomtypes != 13)
	    {
	      fprintf(stderr,"%d is an unallowed number of atoms types ...\nOnly 3, 13 atom types are allowed!\n",atomtypes);
	      exit(1);
	    }
	  i++;
	}
      else if (strcmp(argv[i],"-base")==0)
	{  
	  strcpy(basename,argv[i+1]);
	  strcpy(m[0].filename,basename);
	  i++;
	}
      i++;
    }
  
  //Use the basename with standard extensions.
  if(strcmp(basename,"undef") !=0)
    {
      strcat(m[0].filename,".pdb");
      strcpy(surfacefile,basename);
      strcat(surfacefile,".rsa");
      strcpy(stridefile,basename);
      strcat(stridefile,".stride");
      strcpy(psipredfile,basename);
      strcat(psipredfile,".ss2");
      if(noprof == 0)
	{
	  strcpy(proffile,basename);
	  strcat(proffile,".psi");
	}
    }

  printf("#pdbfile: %s\n",m[0].filename);
  printf("#surface: %s\n",surfacefile);
  printf("#stride:  %s\n",stridefile);
  if(strcmp(proffile,"undef") ==0 || noprof){   
    printf("#profile: %s\n","NOT USED");
  } else {
    printf("#profile: %s\n",proffile);
  }
  printf("#residuetypes: %d\n",residuetypes);
  printf("#atomtypes: %d\n",atomtypes);
  printf("#windowsize: %d\n",windowsize);
  printf("#residuecutoff: %5.3lf\n",sqrt(cutoff));
  printf("#atomcutoff: %5.3lf\n",sqrt(atomcutoff));

  
  //printf("%s\n%s\n%s\n%s\n%s\n%d %d %8.3lf %8.3lf\n",m[0].filename,surfacefile,stridefile,psipredfile,proffile,noprof,windowsize,cutoff,atomcutoff);
  //exit(1);
  
  //for(i=0;i<20;i++)
  //  {
  //    printf("%d %c %d %d\n",i,profile_index_to_aa(i),get_res6_no_pointer(profile_index_to_aa(i)),1);
  //  }
  //exit(1);
  if(strcmp(m[0].filename,"0") == 0 ||
     strcmp(surfacefile,"undef") == 0 ||
     strcmp(stridefile,"undef") == 0 ||
     strcmp(psipredfile,"undef") == 0)
    {
      usage();
    }
  
  if(getenv("PROQRESDIR")==NULL)
    {
      fprintf(stderr,"ENV variable PROQRESDIR needs to be set properly\n");
      exit(1);
    }
  strcpy(net1,getenv("PROQRESDIR"));
  strcpy(net2,net1);
  strcpy(net3,net1);
  strcpy(net4,net1);
  strcpy(net5,net1);
  printf("#$PROQRESDIR=%s\n",net1);
  read_net((char*)strcat(net1,"/net1"),&net[0]);
  read_net((char*)strcat(net2,"/net2"),&net[1]);
  read_net((char*)strcat(net3,"/net3"),&net[2]);
  read_net((char*)strcat(net4,"/net4"),&net[3]);
  read_net((char*)strcat(net5,"/net5"),&net[4]);

  //  exit(1);
  //printf("%s %s %s %s %d\n",m[0].filename,surfacefile,stridefile,psipredfile,windowsize%2);
  //prof[0][0]=1;
  //prof[0][1]=2;
  //prof[1][0]=3;
  //printf("%d\n", prof[0][0]);


  //printf("%d \n", strlen(profseq));
  //exit(1);
  //for(i=0;i<strlen(profseq);i++)
  //  {
  //    printf("%c ",profseq[i]);
  //    for(j=0;j<22;j++)
  //	printf("%3.2lf ",prof[i][j]);
  //    printf("\n");
  //  }
  //exit(1);
  stride=read_stride(stridefile);
  psipred=read_psipred2(psipredfile,coil,helix,sheet);
  //printf("%s",psipred);
  if(windowsize%2!=1)
    {
      fprintf(stderr,"Windowsize should be an odd number!\n");
      exit(1);
    }
  //exit(1);
  fp=fopen(surfacefile,"r");
  if(fp!=NULL)
    {
      i=0;
      while(fgets(buff,1000,fp)!=NULL)
	{
	  sscanf(buff,"%s",line_flag);
	  if(strcmp("RES",line_flag)==0)
	    {
	      strncpy_NULL(temp_surface,&buff[36],5);
	      surfaces[i]=atof(temp_surface);
	      i++;
	    }
	}
      res_in_surffile=i;
      //for(j=0;j<i;j++)
      //{
      //  printf("%lf\n",surfaces[j]);
      //  }
      // exit(1);

    }
  else
    {
      fprintf(stderr,"Error cannot open %s\n",surfacefile);
      exit(1);
    }
  //exit(1);
  //cutoff=strtod(argv[2],(char**)(argv[2]+strlen(argv[2])));
  //cutoff=49; //cutoff*cutoff (7A for residues);
  //atomcutoff=25; //(5A for atom contacts)
  if(read_molecules(m,'a')==0)
    { 
      if(strcmp(proffile,"undef")!=0 && !noprof)
	{
	  printf("#Reading prof...\n");
	  read_profile(proffile,prof,profseq);

	  // Check if profile has all zeros column.
	  //printf("%s\n",profseq);
	  for(i=0;i<strlen(profseq);i++)
	    {
	      check=0;
	      for(j=0;j<20;j++)
		{
		  //  printf("%lf ",prof[i][j]);
		  if(prof[i][j]!=0)
		    check=1;
		}
	      //printf("\n");
	      if(check==0)
		{
		  //  printf("%d %c\n",i+1,profseq[i]);
		  prof[i][aa(profseq[i])]=1;

		}
	    }
	  //printf("%s\n",profseq);
	}
      else //initialize prof and profseq to standard values
	{
	  //	   printf("#Use identity prof...\n");
	  //initialize prof and profseq to standard values
	  strcpy(profseq,m[0].sequence);
	  for(i=0;i<strlen(profseq);i++)
	    {
	      prof[i][aa(profseq[i])]=1;
	    }
	}
      //      exit(1);
 
      for(i=0;i<m[0].atoms;i++)  
	{
	  if(m[0].atm[i].rescount!=current_res_i)
	    {
	      for(j=i,current_res_j=current_res_i;j<m[0].atoms;j++)
		{
		  if(m[0].atm[j].rescount!=current_res_j)
		    {
		      //printf("%d %d %d %d %f %f\n",current_res_j,current_res_i,i,j,crd(m,i,j),crd(m,j,i));
		      //	     {
		      //    printf("%d %d %d %d %f %f\n",m[0].atm[j].rescount,m[0].atm[i].rescount,i,j,crd(m,i,j),crd(m,j,i));
		      // printf("%d %s %s %i\n" ,m[0].atm[i].rescount,m[0].atm[i].name,m[0].atm[i].residue,m[0].atm[i].resnum);
		      //printf("%d %s %s %i\n" ,m[0].atm[j].rescount,m[0].atm[j].name,m[0].atm[j].residue,m[0].atm[j].resnum);
		      //printf("%d %d %d %d %f %f\n",m[0].atm[j].rescount,m[0].atm[i].rescount,i,j,crd(m,i,j),crd(m,j,i));
		      //if(abs(m[0].atm[i].rescount-m[0].atm[j].rescount)>5 &&
		      //if(abs(m[0].atm[i].rescount-m[0].atm[j].rescount)>5 && crd(m,i,j)<cutoff)
		      if(crd(m,i,j)<cutoff)
			{
			  
			  //printf("%d %d %d %d %f %f\n",m[0].atm[j].rescount,m[0].atm[i].rescount,i,j,crd(m,i,j),crd(m,j,i));
			  //printf("%d\n",contacts[current_res_i][current_res_j]);
			  //rescontacts2[current_res_i][current_res_j]=1;
			  //rescontacts2[current_res_j][current_res_i]=1;
			  //rescontacts[current_res_i][get_res6(m[0].atm[i].residue)][get_res6(m[0].atm[j].residue)]++;
			  //rescontacts[current_res_j][get_res6(m[0].atm[j].residue)][get_res6(m[0].atm[i].residue)]++;
			  
			  //res_contact_map[current_res_i][current_res_j]=1;
			  //res_contact_map[current_res_j][current_res_i]=1;
			  
			  if(abs(m[0].atm[i].rescount-m[0].atm[j].rescount)>5)
			    {
			      res_contact_map[current_res_i][current_res_j]=1;
			      res_contact_map[current_res_j][current_res_i]=1;

			      //rescontacts6[current_res_i][get_res6(m[0].atm[j].residue)]++;
			      //rescontacts6[current_res_j][get_res6(m[0].atm[i].residue)]++;
			    }
			  //res_contacts[get_res(m[0].atm[i].residue)][get_res(m[0].atm[j].residue)]++;
			  //tot_res_contacts++;
			}
		      current_res_j++;
		    }
		}				     
	      current_res_i++;
	    }
	}
      for(i=0;i<m[0].atoms;i++)  
	{
	  for(j=i;j<m[0].atoms;j++)  
	    {
	      //distance(m[0],i,j);
	      if(abs(m[0].atm[i].rescount-m[0].atm[j].rescount)>1 && distance(m,i,j)<atomcutoff)
		{
		  current_atomtype_i=get_atomtype(m[0].atm[i].name,m[0].atm[i].residue);
		  current_atomtype_j=get_atomtype(m[0].atm[j].name,m[0].atm[j].residue);

		  current_atomtype3_i=get_atomtype3(m[0].atm[i].name,m[0].atm[i].residue);
		  current_atomtype3_j=get_atomtype3(m[0].atm[j].name,m[0].atm[j].residue);

		  atomcontacts[m[0].atm[i].rescount][current_atomtype_i][current_atomtype_j]++;
		  atomcontacts[m[0].atm[j].rescount][current_atomtype_j][current_atomtype_i]++;

		  if(current_atomtype3_j != 3 || current_atomtype3_i !=3)
		    {
		      atomcontacts3[m[0].atm[i].rescount][current_atomtype3_i][current_atomtype3_j]++;
		      atomcontacts3[m[0].atm[j].rescount][current_atomtype3_j][current_atomtype3_i]++;
		    }
		  //atomcontacts[m[0].atm[i].rescount][current_atomtype_j][current_atomtype_i]++;
		  //atomcontacts[m[0].atm[j].rescount][current_atomtype_i][current_atomtype_j]++;
		  //printf("%d %s %s %d<--> %d %s %s %d=> %lf\n",m[0].atm[i].rescount,m[0].atm[i].name,m[0].atm[i].residue,get_atomtype(m[0].atm[i].name,m[0].atm[i].residue),m[0].atm[j].rescount,m[0].atm[j].name,m[0].atm[j].residue,get_atomtype(m[0].atm[j].name,m[0].atm[j].residue),distance(m,i,j));
		}
	      
	    }
	}
      // exit(1);
      if(res_in_surffile == m[0].residues && strlen(psipred) == strlen(stride) && strlen(stride)==m[0].residues)
	{
	  for(i=-(windowsize-1)/2;i<m[0].residues-(windowsize-1)/2;i++)
	    {
	      // printf("%d %d\n",i,i+windowsize-1);
	      
	      //Over each residue in window 
	      //res=0;
	      total_res_contacts_count=0;
	      total_res_contacts=0;
	      rescount=0;
	      atomcount=0;
	      surface_total=0;
	      surface_total_prof=0;
	      for(j=i;j<i+windowsize;j++)
		{
		  if(j>=0 && j<m[0].residues)
		    {
		      strncpy_NULL(res1,&m[0].sequence[j],1);
		      //printf("%s %c %lf %lf %lf\n",res1,stride[j],surfaces[j],buried_cutoff,exposed_cutoff);
		      
		      if(residuetypes==623)
			restype1=get_residue_type623(res1,stride[j],surfaces[j],buried_cutoff,exposed_cutoff);
		      if(residuetypes==622)
			restype1=get_residue_type622(res1,stride[j],surfaces[j],buried_cutoff);
		      if(residuetypes==632)
			restype1=get_residue_type632(res1,stride[j],surfaces[j],buried_cutoff);
		      if(residuetypes==20)
			restype1=aa(m[0].sequence[j]);
		      if(residuetypes == 6)
			restype1=get_res6(res1);
		      current_res_j=get_res6(res1);
		      //printf("%d %c(%d:%3.2lf)",j,m[0].sequence[j],aa(m[0].sequence[j]),prof[j][aa(m[0].sequence[j])]);
		      // printf("[");
		      //Over contacting residues
		      //Calculate all pair of contacts from the profile vector
		      for(k=0;k<m[0].residues;k++)
			{
			  if(res_contact_map[j][k]==1)
			    {
			      //Contact between j and k rows in profile (between res1 and res2)
			      strncpy_NULL(res2,&m[0].sequence[k],1);
			      if(residuetypes==623)
				restype2=get_residue_type623(res2,stride[k],surfaces[k],buried_cutoff,exposed_cutoff);
			      if(residuetypes==622)
				restype2=get_residue_type622(res2,stride[k],surfaces[k],buried_cutoff);
			      if(residuetypes==632)
				restype2=get_residue_type632(res2,stride[k],surfaces[k],buried_cutoff);
			      if(residuetypes==20)
				restype2=aa(m[0].sequence[k]);
			      if(residuetypes == 6)
				restype2=get_res6(res2);

			      total_res_contacts_count++;

			      //loop over all pairs in profile row j and k
			      if(!noprof)
				{
				  if(residuetypes==20)
				    {
				      for(l=0;l<20;l++) 
					{
					  for(n=0;n<20;n++) 
					    {
					      
					      rescontacts_prof[l][n]=rescontacts_prof[l][n]+prof[j][l]*prof[k][n];
					      total_res_contacts=total_res_contacts+prof[j][l]*prof[k][n];
					    }
					}
				    }
				  else
				    {
				      for(l=0;l<20;l++) 
					{
					  res1[0]=profile_index_to_aa(l);
					  res1[1]='\0';
					  if(residuetypes==623)
					    restype1=get_residue_type623(res1,stride[j],surfaces[j],buried_cutoff,exposed_cutoff);
					  if(residuetypes==622)
					    restype1=get_residue_type622(res1,stride[j],surfaces[j],buried_cutoff);
					  if(residuetypes==632)
					    restype1=get_residue_type632(res1,stride[j],surfaces[j],buried_cutoff);
					  if(residuetypes == 6)
					    restype1=get_res6(res1);
					  for(n=0;n<20;n++) 
					    {
					      res2[0]=profile_index_to_aa(n);
					      res2[1]='\0';
					      if(residuetypes==623)
						restype2=get_residue_type623(res2,stride[k],surfaces[k],buried_cutoff,exposed_cutoff);
					      if(residuetypes==622)
						restype2=get_residue_type622(res2,stride[k],surfaces[k],buried_cutoff);
					      if(residuetypes==632)
						restype2=get_residue_type632(res2,stride[k],surfaces[k],buried_cutoff);
					      if(residuetypes == 6)
						restype1=get_res6(res1);				      
					      // printf("%d %d: %s %s %d %d %d %d\n",i,j,res1,res2,l,n,restype1,restype2);
					      rescontacts_prof[restype1][restype2]=rescontacts_prof[restype1][restype2]+prof[j][l]*prof[k][n];
					      total_res_contacts=total_res_contacts+prof[j][l]*prof[k][n];
					    }
					}
				    }
				}
			      else
				{
				  //printf("No prof %d!\n",noprof);
				  rescontacts_prof[restype1][restype2]+=1.0;
				  //	  printf("No prof %d (%d,%d) %lf!\n",noprof, restype1,restype2,rescontacts_prof[restype1][restype2]);
				  total_res_contacts++;
				}
			      // printf("%d %lf\n",k+1,prof_sum_tmp);
			    }
			}
		      // printf("]");
		      //exit(1);
			  //printf("%c(%d:%3.2lf)",m[0].sequence[j],aa(m[0].sequence[j]),prof[j][aa(m[0].sequence[j])]);
		      //*** RESIDUE CONTACTS ***
		     //for(k=0;k<6;k++)
		     //	{
		     //	  res[current_res_j][k]+=rescontacts[j][k];
		     //	  rescount+=rescontacts[j][k];
		     //	}
		     ////*** ATOM CONTACTS ***
		      for(k=0;k<13;k++)
			{
			  for(l=0;l<13;l++)
			    {
			      atom[k][l]+=atomcontacts[j][k][l];
			      atomcount+=atomcontacts[j][k][l];
			    }
			}
		      for(k=0;k<3;k++)
			{
			  for(l=0;l<3;l++)
			    {
			      atom3[k][l]+=atomcontacts3[j][k][l];
			      atomcount3+=atomcontacts3[j][k][l];
			    }
			}

		      
		      //printf("%lf\n ",surfaces[j]);
		      if(surfaces[j]<=25){
		     	//surface25[current_res_j]++;
		     	surface_total++;
			//printf("SURF25\n");
			for(k=0;k<20;k++)
			  {
			    // printf("%lf ",prof[j][k]);
			    surface25_prof[k]+=prof[j][k];
			    surface_total_prof+=prof[j][k];
			  }
		      }
		      if(surfaces[j]>25 && surfaces[j]<=50){
			//printf("SURF50\n");
		     	//surface50[current_res_j]++;
		     	surface_total++;
			for(k=0;k<20;k++)
			  {
			    surface50_prof[k]+=prof[j][k];
			    surface_total_prof+=prof[j][k];
			  }
		      }
		      if(surfaces[j]>50 && surfaces[j]<=75){
			//printf("SURF75\n");
		     	//surface75[current_res_j]++;
		     	surface_total++;
			for(k=0;k<20;k++)
			  {
			    surface75_prof[k]+=prof[j][k];
			    surface_total_prof+=prof[j][k];
			  }
		      }
		      if(surfaces[j]>75){
			//printf("SURF100\n");
		     	//surface100[current_res_j]++;
		     	surface_total++;
			for(k=0;k<20;k++)
			  {
			    //  printf("%lf ",prof[j][k]);
			    surface100_prof[k]+=prof[j][k];
			    surface_total_prof+=prof[j][k];
			  }
			//printf("\n");

		      }
		      ///printf("SURF: %d %c %lf %d\n",j,m[0].sequence[j],surfaces[j],surface_total);
		      //for(k=0;k<20;k++)
		      //	  {
		      //	    printf("%c %d (%lf %lf %lf %lf)\n",profile_index_to_aa(k),k,surface25_prof[k],surface50_prof[k],surface75_prof[k],surface100_prof[k]);
		      //	  }
		      //printf("\n");
		    }
		  //printf("END\n");
		 
		  //printf("%lf %lf %lf %lf\n",surface25[current_res_j],surface50[current_res_j],surface75[current_res_j],surface100[current_res_j]);
		  //exit(1);
		}
	      //exit(1);
	      //  printf("%d %lf %d %lf\n",total_res_contacts_count,total_res_contacts,surface_total,surface_total_prof);

	      //temp=0;
	      //for(j=0;j<strlen(m[0].sequence);j++)
	      //	{
	      //	  printf("%d: ",j);
	      //	  for(k=0;k<20;k++)
	      //	    {
	      //	      printf("%lf ",prof[j][k]);
	      //	    }
	      //	  printf("\n");
	      //	}
	      //exit(1);
	      //*** RESIDUE CONTACTS ***
	      if(residuetypes==6){
		max_j=6;
		max_k=6;
	      }
	      if(residuetypes==20){
		max_j=20;
		max_k=20;
	      }
	      if(residuetypes==623 || residuetypes==632){
		max_j=36;
		max_k=36;
	      }
	      if(residuetypes==622){
		max_j=24;
		max_k=24;
	      }
	      for(j=0,n=0;j<max_j;j++){
		for(k=j;k<max_k;k++){
		  if(k>=j && total_res_contacts!=0){
		    if(k!=j){
		      //printf("j:%d k:%d %lf %lf\n",j,k,rescontacts_prof[k][j]+rescontacts_prof[j][k],total_res_contacts);
		    resfrac[n]=(rescontacts_prof[k][j]+rescontacts_prof[j][k])/total_res_contacts;
		    }else{
		      // printf("j:%d k:%d %lf %lf\n",j,k,rescontacts_prof[k][j]+rescontacts_prof[j][k],total_res_contacts);
		      resfrac[n]=rescontacts_prof[j][k]/total_res_contacts;
		    }
		  }else{
		    resfrac[n]=0.0;
		  }
		  n++;
		}
	      }
	      
	      	      
	      //printf("%d\n",n);
	      //*** SURFACES *** RESIDUE CONTACTS ***
	      //Keep these for the time being since it sounds strange to use other the other types that are defined using surfaces....
	      tmp_sum=0;
	      if(residuetypes==6)
		{
		  for(j=0;j<20;j++)
		    {
		      res1[0]=profile_index_to_aa(j);
		      res1[1]='\0';
		      current_res_j=get_res6(res1);
		      surface25[current_res_j]+=surface25_prof[j];
		      surface50[current_res_j]+=surface50_prof[j];
		      surface75[current_res_j]+=surface75_prof[j];
		      surface100[current_res_j]+=surface100_prof[j];
		      //  printf("%s %d->%d (%lf %lf %lf %lf)\n",res1,j,current_res_j,surface25_prof[j],surface50_prof[j],surface75_prof[j],surface100_prof[j]);
		      
		    }
		  for(j=0;j<6;j++)
		    {
		      
		      surface25_prof[j]=surface25[j]/surface_total;
		      surface50_prof[j]=surface50[j]/surface_total;
		      surface75_prof[j]=surface75[j]/surface_total;
		      surface100_prof[j]=surface100[j]/surface_total;
		      
		    }

		}
	      else
		{
		  for(j=0;j<20;j++)
		    {
		      surface25_prof[j]/=surface_total;
		      surface50_prof[j]/=surface_total;
		      surface75_prof[j]/=surface_total;
		      surface100_prof[j]/=surface_total;
		    }
		}

	      
	      //	      printf("Sum-surf: %lf %d %lf\n",tmp_sum,surface_total,surface_total_prof);
	      // for(j=0;j<6;j++)
	      //	{
	      //	  printf("%lf ", surface25_prof[j]);
	      //	}
	      //printf("\n";)

		//exit(1);
	      //for(j=0;j<6;j++)
	      //	{
	      //	  surface25[j]=surface25[j]/surface_total; //m[0].residues;
	      //	  surface50[j]=surface50[j]/surface_total; //m[0].residues;
	      //	  surface75[j]=surface75[j]/surface_total; //m[0].residues;
	      //	  surface100[j]=surface100[j]/surface_total; //m[0].residues;
	      //	}
	      //END OF ONE WINDOW

	      //*** RESIDUE CONTACTS ***
	     //for(k=0,n=0;k<6;k++){
	     //	//if(total_res_contacts!=0){
	     //	//  resfrac[n]=(double)res[k][k]/rescount;
	     //	//}else{
	     //	//  resfrac[n]=0.0;
	     //	//}
	     //	//n++;
	     //	//for(l=k+1;l<6;l++){
	     //	//  if(rescount!=0){
	     //	//    resfrac[n]=(double)(res[k][l]+res[l][k])/rescount;
	     //	//    //printf("%d\t",res[k][l]+res[l][k]);
	     //	//  }else{
	     //	//    resfrac[n]=0;
	     //	//  }
	     //	//  n++;
	     //	//}
	     //	//printf("%d\n",surface_total);
	     //	surface25[k]=surface25[k]/surface_total; //m[0].residues;
	     //	surface50[k]=surface50[k]/surface_total; //m[0].residues;
	     //	surface75[k]=surface75[k]/surface_total; //m[0].residues;
	     //	surface100[k]=surface100[k]/surface_total; //m[0].residues;
	     //
	     //}

	      //for(k=0,n=0;k<13;k++){
	      //	for(l=0;l<13;l++){
	      //	  //printf("%d ",atom[k][l]);
	      //	}
	      //	//printf("\n");
	      //}
	      //printf("%d\n",atomcount);
		
	      //*** ATOM CONTACTS ***
	      for(k=0,n=0;k<13;k++){
		if(atomcount!=0){
		  atomfrac[n]=(double)atom[k][k]/atomcount;
		}else{
		  atomfrac[n]=0;
		}
		n++;
		for(l=k+1;l<13;l++){
		  if(atomcount!=0){
		    atomfrac[n]=(double)(atom[k][l]+atom[l][k])/atomcount;
		    //printf("%d\t",atom[k][l]+atom[l][k]);
		  }else{
		    atomfrac[n]=0;
		  }
		  n++;
		}
	      }

	      for(k=0,n=0;k<3;k++){
		if(atomcount3!=0){
		  atomfrac3[n]=(double)atom3[k][k]/atomcount3;
		}else{
		  atomfrac3[n]=0;
		}
		n++;
		for(l=k+1;l<3;l++){
		  if(atomcount3!=0){
		    atomfrac3[n]=(double)(atom3[k][l]+atom3[l][k])/atomcount3;
		    //printf("%d\t",atom[k][l]+atom[l][k]);
		  }else{
		    atomfrac3[n]=0;
		  }
		  n++;
		}
	      }
	      
	      //**** SECONDARY STRUCTURE ****
	      current_res_i=i+(windowsize-1)/2;
	      ss=0;
	      if(stride[current_res_i] == 'H')
		{
		  ss=helix[current_res_i];
		}else{
		  if(stride[current_res_i] == 'E')
		    {
		      ss=sheet[current_res_i];
		    }
		  else
		    {
		      ss=coil[current_res_i];
		    }
		}
	    
	      
	      // fprintf(stderr,"%d %c %lf %c %c\n",current_res_i,m[0].sequence[current_res_i],ss,psipred[current_res_i],stride[current_res_i]);
	 

	  	  	      
	      //Put all in the parameters vector
	      l=0;
	      for(k=0;k<91;k++,l++)
		parameters[l]=atomfrac[k];
	      for(k=0;k<21;k++,l++)
		parameters[l]=resfrac[k];
	      for(k=0;k<6;k++,l++)
		parameters[l]=surface25_prof[k];
	      for(k=0;k<6;k++,l++)
		parameters[l]=surface50_prof[k];
	      for(k=0;k<6;k++,l++)
		parameters[l]=surface75_prof[k];
	      for(k=0;k<6;k++,l++)
		parameters[l]=surface100_prof[k];
	      parameters[l]=ss;
	      predicted_quality=0;
	      for(k=0;k<5;k++){
		quality[k]=netfwd(parameters,&net[k]);
		predicted_quality+=quality[k];
	      }
	      predicted_quality/=5;
		
		//      predicted_quality=netfwd(parameters
	      //printf("A=[");
	      //for(k=0;k<137;k++)
	      //	printf("%lf ",parameters[k]);
	      //printf("]\n");
	      printf("%6d %c %8.5lf\n",i+(windowsize-1)/2+1,m[0].sequence[i+(windowsize-1)/2],predicted_quality);
	      //printf("%6d %c %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf\n",i+(windowsize-1)/2+1,m[0].sequence[i+(windowsize-1)/2],predicted_quality, quality[0],quality[1],quality[2],quality[3],quality[4]);
	
	      //Put things to zero again, perhaps an easier way.
	       for(k=0;k<36;k++){
		 for(l=0;l<36;l++){
		   rescontacts_prof[k][l]=0;
		 }
	       }
	       for(k=0;k<20;k++){
		 surface25_prof[k]=0;
		 surface50_prof[k]=0;
		 surface75_prof[k]=0;
		 surface100_prof[k]=0;
	       }
	       for(k=0;k<666;k++)
		 resfrac[k]=0;
	       for(k=0;k<91;k++)
		 atomfrac[k]=0;
	       
	      for(k=0;k<6;k++){
		for(l=0;l<6;l++){
		  res[k][l]=0;
		}
		surface25[k]=0;
		surface50[k]=0;
		surface75[k]=0;
		surface100[k]=0;
	      }
	      for(k=0;k<14;k++){
		for(l=0;l<14;l++){
		  atom[k][l]=0;
		}
	      }
	      
	    }
	  
	  //exit(1);

	}
	  else
	{
	  fprintf(stderr,"Different number of residues in one of these filename (#number of residues)\nsurface file: %s (%d)\npdbfile: %s (%d)\nstridefile: %s (%d)\npsipredfile: %s (%d)\n",surfacefile,res_in_surffile,m[0].filename,m[0].residues,stridefile,strlen(stride),psipredfile,strlen(psipred));
	  printf("%s\n",psipred);
	}

      //     for(i=0;i<m[0].residues;i++)  
      //	{
      //	  printf("%d %c\n",i,m[0].sequence[i]);
      //	  tmp=0;
      //	  //printf("C\tN\tO\tCA\tCH3\tCH/CH2\tC(OO-)\tNH\tNH2\t(C)00-\t=0\tOH\tS\tOXT\n");
      //	  for(j=0;j<6;j++)  
      //	    {
      //	      for(k=0;k<6;k++)  
      //		{
      //		  printf("%d\t",rescontacts[i][j][k]);
      //		}
      //	      printf("\n");
      //	    }
      //	  printf("\n\n");
      //	}
      //     exit(1);

    }
  free(stride);
}


int usage()
{  
  fprintf(stderr,"\nUsage:\nProQres\t args\n\t\t -base [basename (*.ss2,*.rsa,*.stride,*.pdb)]\n\t\t -pdb  [pdb file]\n\t\t -surf [surface file]\n\t\t -stride  [stride file]\n\t\t -psipred [psipred file]\n\t\t\n");
  //fprintf(stderr,"\nUsage:\nProQres\t args\n\t\t -base [basename]\n\t\t -pdb  [pdb file]\n\t\t -surf [surface file]\n\t\t -prof [psiblast profile]\n\t\t -stride  [stride file]\n\t\t -psipred [psipred file]\n\t\t -rc [residue cutoff]\n\t\t -ac [atom cutoff]\n\t\t -bc [buried cutoff (used strict (above is exposed) for 632 and 622)]\n\t\t -ec [exposed cutoff (only 623)]\n\t\t -w  [window size]\n\t\t -t  [residue types 6/20/623/632/622 (623 = 6 res + 2 ss + 3 surfaces)]\n"); //\t\t -w [window size]\n\t\t -rc [residue cutoff (7Å)]\n\t\t -ac [atom cutoff (5Å)]\n");//\t\t -ss2 [calculate ss in another way]\n");
  exit(1);
}     



     
int get_residue_type623(char* res,char ss,double surface,double buried_cutoff,double exposed_cutoff)
{


  int residue_class=0; // The ordinary 6 classes: 0:+,1:-,2:aro,3:polar,4:hydrophobic,:5:special
  int ss_class=0; //0 loop, 1 non-loop
  int surface_class=0; //0 <buried_cutoff, 1 buried_cutoff-exposed_cutoff, 2 > exposed_cutoff
  int index=0;

  if(strcmp("ARG",res)==0 || strcmp("LYS",res)==0)
    residue_class = 0;
  if(strcmp("ASP",res)==0 || strcmp("GLU",res)==0)
    residue_class = 1;
  if(strcmp("HIS",res)==0 || strcmp("PHE",res)==0 || strcmp("TRP",res)==0 || strcmp("TYR",res)==0)
    residue_class = 2;
  if(strcmp("ASN",res)==0 || strcmp("GLN",res)==0 || strcmp("SER",res)==0 || strcmp("THR",res)==0)
    residue_class = 3;
  if(strcmp("ALA",res)==0 || strcmp("ILE",res)==0 || strcmp("LEU",res)==0 || strcmp("MET",res)==0 || strcmp("VAL",res)==0 || strcmp("CYS",res)==0)
    residue_class = 4;
  if(strcmp("GLY",res)==0 || strcmp("PRO",res)==0)
    residue_class = 5;
  if(strcmp("R",res)==0 || strcmp("K",res)==0)
    residue_class = 0;
  if(strcmp("D",res)==0 || strcmp("E",res)==0)
    residue_class = 1;
  if(strcmp("H",res)==0 || strcmp("F",res)==0 || strcmp("W",res)==0 || strcmp("Y",res)==0)
    residue_class = 2;
  if(strcmp("N",res)==0 || strcmp("Q",res)==0 || strcmp("S",res)==0 || strcmp("T",res)==0)
    residue_class = 3;
  if(strcmp("A",res)==0 || strcmp("I",res)==0 || strcmp("L",res)==0 || strcmp("M",res)==0 || strcmp("V",res)==0 || strcmp("C",res)==0)
    residue_class = 4;
  if(strcmp("G",res)==0 || strcmp("P",res)==0)
    residue_class = 5;
  
  if(ss == 'H' || ss == 'E')
    ss_class=1;
  else
    ss_class=0;
  
  if(surface < buried_cutoff)
    surface_class=0;
  if(surface >= buried_cutoff && surface < exposed_cutoff)
     surface_class=1;
  if(surface > exposed_cutoff)
     surface_class=2;


  if(ss_class==0 && surface_class==0)
    index=residue_class;
  if(ss_class==0 && surface_class==1)
    index=residue_class+6;
  if(ss_class==0 && surface_class==2)
    index=residue_class+6*2;
  if(ss_class==1 && surface_class==0)
    index=residue_class+6*3;
  if(ss_class==1 && surface_class==1)
    index=residue_class+6*4;
  if(ss_class==1 && surface_class==2)
    index=residue_class+6*5;
  
  //printf("623: %s %c %lf (%3.2lf %3.2lf) -> %2d %2d %2d -> %d\n",res,ss,surface,buried_cutoff,exposed_cutoff,residue_class,ss_class,surface_class,index); 
  return index;

}
int get_residue_type622(char* res,char ss,double surface,double buried_cutoff)
{


  int residue_class=0; // The ordinary 6 classes: 0:+,1:-,2:aro,3:polar,4:hydrophobic,:5:special
  int ss_class=0; //0 loop, 1 non-loop
  int surface_class=0; //0 <buried_cutoff, 1 buried_cutoff-exposed_cutoff, 2 > exposed_cutoff
  int index=0;

  if(strcmp("ARG",res)==0 || strcmp("LYS",res)==0)
    residue_class = 0;
  if(strcmp("ASP",res)==0 || strcmp("GLU",res)==0)
    residue_class = 1;
  if(strcmp("HIS",res)==0 || strcmp("PHE",res)==0 || strcmp("TRP",res)==0 || strcmp("TYR",res)==0)
    residue_class = 2;
  if(strcmp("ASN",res)==0 || strcmp("GLN",res)==0 || strcmp("SER",res)==0 || strcmp("THR",res)==0)
    residue_class = 3;
  if(strcmp("ALA",res)==0 || strcmp("ILE",res)==0 || strcmp("LEU",res)==0 || strcmp("MET",res)==0 || strcmp("VAL",res)==0 || strcmp("CYS",res)==0)
    residue_class = 4;
  if(strcmp("GLY",res)==0 || strcmp("PRO",res)==0)
    residue_class = 5;
  if(strcmp("R",res)==0 || strcmp("K",res)==0)
    residue_class = 0;
  if(strcmp("D",res)==0 || strcmp("E",res)==0)
    residue_class = 1;
  if(strcmp("H",res)==0 || strcmp("F",res)==0 || strcmp("W",res)==0 || strcmp("Y",res)==0)
    residue_class = 2;
  if(strcmp("N",res)==0 || strcmp("Q",res)==0 || strcmp("S",res)==0 || strcmp("T",res)==0)
    residue_class = 3;
  if(strcmp("A",res)==0 || strcmp("I",res)==0 || strcmp("L",res)==0 || strcmp("M",res)==0 || strcmp("V",res)==0 || strcmp("C",res)==0)
    residue_class = 4;
  if(strcmp("G",res)==0 || strcmp("P",res)==0)
    residue_class = 5;
  
  if(ss == 'H' || ss == 'E')
    ss_class=1;
  else
    ss_class=0;
  
  if(surface < buried_cutoff)
    surface_class=0;
  if(surface >= buried_cutoff)
     surface_class=1;
 
  if(ss_class==0 && surface_class==0)
    index=residue_class;
  if(ss_class==0 && surface_class==1)
    index=residue_class+6;
  if(ss_class==1 && surface_class==0)
    index=residue_class+6*2;
  if(ss_class==1 && surface_class==1)
    index=residue_class+6*3;
  
  //printf("622: %s %c %lf (%3.2lf) -> %2d %2d %2d -> %d\n",res,ss,surface,buried_cutoff,residue_class,ss_class,surface_class,index); 
  return index;

}
int get_residue_type632(char* res,char ss,double surface,double buried_cutoff)
{


  int residue_class=0; // The ordinary 6 classes: 0:+,1:-,2:aro,3:polar,4:hydrophobic,:5:special
  int ss_class=0; //0 loop, 1 helix,2 sheet
  int surface_class=0; //0 <buried_cutoff,  1 > exposed_cutoff
  int index=0;

  if(strcmp("ARG",res)==0 || strcmp("LYS",res)==0)
    residue_class = 0;
  if(strcmp("ASP",res)==0 || strcmp("GLU",res)==0)
    residue_class = 1;
  if(strcmp("HIS",res)==0 || strcmp("PHE",res)==0 || strcmp("TRP",res)==0 || strcmp("TYR",res)==0)
    residue_class = 2;
  if(strcmp("ASN",res)==0 || strcmp("GLN",res)==0 || strcmp("SER",res)==0 || strcmp("THR",res)==0)
    residue_class = 3;
  if(strcmp("ALA",res)==0 || strcmp("ILE",res)==0 || strcmp("LEU",res)==0 || strcmp("MET",res)==0 || strcmp("VAL",res)==0 || strcmp("CYS",res)==0)
    residue_class = 4;
  if(strcmp("GLY",res)==0 || strcmp("PRO",res)==0)
    residue_class = 5;
  if(strcmp("R",res)==0 || strcmp("K",res)==0)
    residue_class = 0;
  if(strcmp("D",res)==0 || strcmp("E",res)==0)
    residue_class = 1;
  if(strcmp("H",res)==0 || strcmp("F",res)==0 || strcmp("W",res)==0 || strcmp("Y",res)==0)
    residue_class = 2;
  if(strcmp("N",res)==0 || strcmp("Q",res)==0 || strcmp("S",res)==0 || strcmp("T",res)==0)
    residue_class = 3;
  if(strcmp("A",res)==0 || strcmp("I",res)==0 || strcmp("L",res)==0 || strcmp("M",res)==0 || strcmp("V",res)==0 || strcmp("C",res)==0)
    residue_class = 4;
  if(strcmp("G",res)==0 || strcmp("P",res)==0)
    residue_class = 5;
  
  
  if(ss == 'H') 
    ss_class=1;
  if(ss == 'E')
    ss_class=2;
  //OTHERWISE ss_class=0; by default.
  
  if(surface < buried_cutoff)
    surface_class=0;
  if(surface >= buried_cutoff)
     surface_class=1;
 
  if(ss_class==0 && surface_class==0)
    index=residue_class;
  if(ss_class==0 && surface_class==1)
    index=residue_class+6;
  if(ss_class==1 && surface_class==0)
    index=residue_class+6*2;
  if(ss_class==1 && surface_class==1)
    index=residue_class+6*3;
  if(ss_class==2 && surface_class==0)
    index=residue_class+6*4;
  if(ss_class==2 && surface_class==1)
    index=residue_class+6*5;
  //printf("632: %s %c %lf (%3.2lf) -> %2d %2d %2d -> %d\n",res,ss,surface,buried_cutoff,residue_class,ss_class,surface_class,index); 
  return index;

}
