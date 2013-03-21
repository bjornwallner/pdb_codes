#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.h"

main(argc,argv)		/* Main routine */
     int argc;
     char *argv[];
{
  molecule      m[1];
  int           i,j,k,l,n;
  double        dist;
  double        mean=0;
  int           atoms=0;
  int           residues=0;
  int           rescontacts[MAXRES][6]={{}};
  int           atomcontacts[MAXRES][14][14]={{{}}};
  int           atomtotal=0;
  int           restotal=0;
  int           current_res_i=0;
  int           current_res_j=0;
  int           current_atomtype_i=0;
  int           current_atomtype_j=0;
  int           alternative_ss_calc=0;
  int           first=0;
  double        cutoff=49; //residue cutoff 7
  double        atomcutoff=25; // atom cutoff 5;
  double        surfaces[MAXRES]={};
  int           surface_total=0;
  double        surface25[6]={};
  double        surface50[6]={};
  double        surface75[6]={};
  double        surface100[6]={};
  double     resfrac[21]={};
  double     atomfrac[91]={};
  double     ss=0;
  int        res[6][6]={{}};
  int        atom[14][14]={{}};
  int        rescount=0;
  int        atomcount=0;
  int        windowsize;
  FILE          *fp;
  char          surfacefile[500]="undef";
  char          stridefile[500]="undef";
  char          psipredfile[500]="undef";
  double        helix[MAXRES];
  double        coil[MAXRES];
  double        sheet[MAXRES];
  int           res_in_surffile=0;
  char          buff[1000];
  char          line_flag[10];
  char          temp_surface[10];
  char          temp_res[2];
  char *stride;
  char *psipred;
  double temp;
  int tmp=0;
  //float        sum=0;
  /* Parse command line for PDB filename */
  i=1;

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
//else if (strcmp(argv[i],"-ss2")==0)
//	{  
//	  alternative_ss_calc=1;
//	}
      i++;
    }
  //  printf("%s %s %s %s %d %8.3lf %8.3lf\n",m[0].filename,surfacefile,stridefile,psipredfile,windowsize,cutoff,atomcutoff);
  
  if(strcmp(m[0].filename,"0") == 0 ||
     strcmp(surfacefile,"undef") == 0 ||
     strcmp(stridefile,"undef") == 0 ||
     strcmp(psipredfile,"undef") == 0)
    {
      usage();
    }
  
  //printf("%s %s %s %s %d\n",m[0].filename,surfacefile,stridefile,psipredfile,windowsize%2);
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
		      if(abs(m[0].atm[i].rescount-m[0].atm[j].rescount)>5 && crd(m,i,j)<cutoff)
			{
			  //printf("%d %d %d %d %f %f\n",m[0].atm[j].rescount,m[0].atm[i].rescount,i,j,crd(m,i,j),crd(m,j,i));
			  //printf("%d\n",contacts[current_res_i][current_res_j]);
			  //rescontacts2[current_res_i][current_res_j]=1;
			  //rescontacts2[current_res_j][current_res_i]=1;
			  //rescontacts[current_res_i][get_res6(m[0].atm[i].residue)][get_res6(m[0].atm[j].residue)]++;
			  //rescontacts[current_res_j][get_res6(m[0].atm[j].residue)][get_res6(m[0].atm[i].residue)]++;
			  rescontacts[current_res_i][get_res6(m[0].atm[j].residue)]++;
			  rescontacts[current_res_j][get_res6(m[0].atm[i].residue)]++;
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

		  atomcontacts[m[0].atm[i].rescount][current_atomtype_i][current_atomtype_j]++;
		  atomcontacts[m[0].atm[j].rescount][current_atomtype_j][current_atomtype_i]++;
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
	      //Over each residue in window 
	      //res=0;
	      rescount=0;
	      atomcount=0;
	      surface_total=0;
	      for(j=i;j<i+windowsize;j++)
		{
		  if(j>=0 && j<m[0].residues)
		    {
		      strncpy_NULL(temp_res,&m[0].sequence[j],1);
		      current_res_j=get_res6(temp_res);
		      //printf("%d %c %d\n",j,m[0].sequence[j],get_res6(temp_res));
		      //*** RESIDUE CONTACTS ***
		      for(k=0;k<6;k++)
			{
			  res[current_res_j][k]+=rescontacts[j][k];
			  rescount+=rescontacts[j][k];
			}
		      //*** ATOM CONTACTS ***
		      for(k=0;k<13;k++)
			{
			  for(l=0;l<13;l++)
			    {
			      atom[k][l]+=atomcontacts[j][k][l];
			      atomcount+=atomcontacts[j][k][l];
			    }
			}
		      //printf("%lf ",surfaces[j]);
		      if(surfaces[j]<25){
			surface25[current_res_j]++;
			surface_total++;
		      }
		      if(surfaces[j]>25 && surfaces[j]<=50){
			surface50[current_res_j]++;
			surface_total++;
		      }
		      if(surfaces[j]>50 && surfaces[j]<=75){
			surface75[current_res_j]++;
			surface_total++;
		      }
		      if(surfaces[j]>75){
			surface100[current_res_j]++;
			surface_total++;
		      }
		    }
		  //printf("%lf %lf %lf %lf\n",surface25[current_res_j],surface50[current_res_j],surface75[current_res_j],surface100[current_res_j]);
		}
	      //printf("\n");
	      

	      //END OF ONE WINDOW

	      //*** RESIDUE CONTACTS ***
	      for(k=0,n=0;k<6;k++){
		if(rescount!=0){
		  resfrac[n]=(double)res[k][k]/rescount;
		}else{
		  resfrac[n]=0;
		}
		n++;
		for(l=k+1;l<6;l++){
		  if(rescount!=0){
		    resfrac[n]=(double)(res[k][l]+res[l][k])/rescount;
		    //printf("%d\t",res[k][l]+res[l][k]);
		  }else{
		    resfrac[n]=0;
		  }
		  n++;
		}
		//printf("%d\n",surface_total);
		surface25[k]=surface25[k]/surface_total; //m[0].residues;
		surface50[k]=surface50[k]/surface_total; //m[0].residues;
		surface75[k]=surface75[k]/surface_total; //m[0].residues;
		surface100[k]=surface100[k]/surface_total; //m[0].residues;

	      }

	      for(k=0,n=0;k<13;k++){
		for(l=0;l<13;l++){
		  //printf("%d ",atom[k][l]);
		}
		//printf("\n");
	      }
	      //printf("%d\n",atomcount);
		
	      //*** ATOM CONTACTS ***
	      for(k=0,n=0;k<13;k++){
		if(atomcount!=0){
		  atomfrac[n]=(double)atom[k][k]/atomcount;
		}else{
		  resfrac[n]=0;
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


	      
	      for(k=0;k<91;k++)
		printf("%lf ",atomfrac[k]);
	      for(k=0;k<21;k++)
		printf("%lf ",resfrac[k]);
	      for(k=0;k<6;k++)
		printf("%lf ",surface25[k]);
	      for(k=0;k<6;k++)
		printf("%lf ",surface50[k]);
	      for(k=0;k<6;k++)
		printf("%lf ",surface75[k]);
	      for(k=0;k<6;k++)
		printf("%lf ",surface100[k]);
	      printf("%lf\n",ss);
	      // printf("\n\n=================\nNext window\n=================\n\n");
	      
	      

	      //Put things to zero again, perhaps an easier way.
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
	  


	}
      else
	{
	  fprintf(stderr,"Different number of residues in one of these filename (#number of residues)\nsurface file: %s (%d)\npdbfile: %s (%d)\nstridefile: %s (%d)\npsipredfile: %s (%d)\n",surfacefile,res_in_surffile,m[0].filename,m[0].residues,stridefile,strlen(stride),psipredfile,strlen(psipred));
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
  fprintf(stderr,"\nUsage:\ncalc_input2\t args\n\t\t -pdb [pdb file]\n\t\t -surf [surface file]\n\t\t -stride [stride file]\n\t\t -psipred [psipred file]\n\t\t -w [window size]\n\t\t -rc [residue cutoff (7Å)]\n\t\t -ac [atom cutoff (5Å)]\n");//\t\t -ss2 [calculate ss in another way]\n");
  exit(1);
}     

     
