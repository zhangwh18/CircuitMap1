#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#define Pi 3.1415926535897932384626433832795

#define Eexp -20.0
#define L 8
#define M 8
#define N 128

#define DN 60
#define checknum 10000
#define ensemble 2
double zt[checknum+1][N],zvt[checknum+1][N],HFt[checknum+1][N];
double complex Pkt[checknum+1][N],Qkt[checknum+1][N],V[N][N];
double H_bar[checknum+1],logHFt[checknum+1][N],Sigma1[ensemble][checknum+1]={0},Xi[ensemble][checknum+1]={0};
double logEktbar[checknum+1][ensemble],Hkbar[checknum+1],Sigma2[ensemble][checknum+1]={0};


double potentials(double Dz[],double parameter[])
{
	double U;
	U=parameter[2]*(Dz[0]+Dz[1]+Dz[2])*(Dz[0]+Dz[1]+Dz[2])\
		+0.125*parameter[0]*(Dz[0]*Dz[0]*Dz[0]+Dz[1]*Dz[1]*Dz[1]*Dz[1]+Dz[2]*Dz[2]*Dz[2]*Dz[2])\
		+parameter[1]*(Dz[0]*Dz[1]*Dz[0]*Dz[1]+Dz[1]*Dz[2]*Dz[1]*Dz[2]+Dz[0]*Dz[2]*Dz[0]*Dz[2]);
	return U;
}


double Forces(double Dz[],double parameter[])
{
	double forcesT;
	forcesT=-2.*parameter[2]*(4*(Dz[0]+Dz[1]+Dz[2])+Dz[3]+Dz[4]+Dz[5]+Dz[6]+Dz[7]+Dz[8])\
			-0.5*parameter[0]*(Dz[0]*Dz[0]*Dz[0]+Dz[1]*Dz[1]*Dz[1]+Dz[2]*Dz[2]*Dz[2])\
			-2.*parameter[1]*(Dz[0]*Dz[0]*(Dz[1]+Dz[2])+Dz[1]*Dz[1]*(Dz[2]+Dz[0])+Dz[2]*Dz[2]*(Dz[0]+Dz[1])\
			+Dz[0]*(Dz[3]*Dz[3]+Dz[4]*Dz[4])+Dz[1]*(Dz[5]*Dz[5]+Dz[6]*Dz[6])+Dz[2]*(Dz[7]*Dz[7]+Dz[8]*Dz[8]));
	return forcesT;
}

int main()
{		
	clock_t start_run,finish_run;	
	start_run = clock();
	double time_count;
	double a0=1.421E-10,alpha0=155.9,beta0=25.5,gamma0=7.4,m0=1.9926465384E-26,t0=5E-16;
	double t,m,a,alpha,beta,gamma;
	double a1[2],a2[2],b1[2],b2[2],R[3][2],parameter[3];
	double epsinlon,Etot,Erate,Escaled,Epsinlon;
	int i,j,l,n,ii,jj,kk,ll,mm,rr,N_excite;
	double x[N],y[N],r[N][2];
	int nn[N][3],nnn[N][6];
	double kslop=0,kb=0,k_x1[M][L],k_y1[M][L],k[N][2],koder[1][2],f_k[N/2],omega2_1[N/2],omega2_2[N/2];
	double omega_1[N/2],omega_2[N/2],omega_1oder,omega_2oder,omega[N],Omega2[N]; 
	double complex e[2][2][N/2],pp1[N],pp2[N];
	double den[N/2],num[N/2],phi[N/2];
	double qt0[N],pt0[N],qt1[N],qt2[N],vt[N],pt[N];
	double Pk0[N],Qk0[N]; 
	double complex q0[N],p0[N],q[N],p[N]; 
	double qi,qj1,qj2,qj3,qk1,qk2,qk3,qk4,qk5,qk6;
	double P_total=0,epsilon_real,lambda;
	double rand_E[N_excite],pE[N_excite],Eg[N],phi_r[N],H_total;
	double force[N],T[N],U[N],Dz[9]={0};
	double SSigma1,SSigma2,sumEktbar,kstart,Ektbar[N],eta,Wk[N];
	
	int num_of_steps=0;
	int min_excited_mode,max_excited_mode;
	int iteration,N_excited_modes;
	
	//
	double TimeT[checknum];
	//
	epsinlon=pow(10,Eexp);
	Etot=N*epsinlon;
	Erate=pow(t0,2)/(m0*pow(a0,2));
	Escaled=Erate*Etot;Epsinlon=Escaled/N;
//	N_excited_modes=floor(0.1*N);
	min_excited_mode=2-1;
//max_excited_mode=min_excited_mode+N_excited_modes-1;
	max_excited_mode=19-1;
	N_excited_modes=max_excited_mode-min_excited_mode+1;
//	printf("%d %d %d\n",N_excited_modes,min_excited_mode,max_excited_mode);
	printf("Nunber_of excited_modes:%d\n",N_excited_modes);

	printf("interval_t0=%.6e s\n",t0);
	t=1;m=1;a=1;
	alpha=pow(t0,2)/m0*alpha0;beta=pow(t0,2)/m0*beta0;gamma=pow(t0,2)/m0*gamma0;
	
	parameter[0]=alpha;parameter[1]=beta;parameter[2]=gamma;
	
	a1[0]=3*a/2;a1[1]=sqrt(3)*a/2;a2[0]=3*a/2;a2[1]=-sqrt(3)*a/2;
	b1[0]=2*Pi/(3*a);b1[1]=sqrt(3)*2*Pi/(3*a);b2[0]=2*Pi/(3*a);b2[1]=-sqrt(3)*2*Pi/(3*a);
	R[0][0]=(a1[0]+a2[0])/3;R[0][1]=(a1[1]+a2[1])/3;
	R[1][0]=(a1[0]-2*a2[0])/3;R[1][1]=(a1[1]-2*a2[1])/3;
	R[2][0]=(a2[0]-2*a1[0])/3;R[2][1]=(a2[1]-2*a1[1])/3;
	
	for(jj=1;jj<=L;jj++)
	{
		for(ii=1;ii<=M;ii++)
		{
			for(n=1;n<=2;n++)
			{
				j=2*((jj-1)*M+ii-1)+n;	
				x[j-1]=(jj-1)*a1[0]+(ii-1)*a2[0]+(n)*a;
				y[j-1]=(jj-1)*a1[1]+(ii-1)*a2[1];
				r[j-1][0]=x[j-1];
				r[j-1][1]=y[j-1];
				/*printf("%d %d %d\n",n,jj,ii);*/
				if(n==1) 
				{
					nn[j-1][0]=j+1;	
					if(jj==1) nn[j-1][1]=j+2*(L-1)*M+1;
					else 	  nn[j-1][1]=j-2*M+1;
					if(ii==1) nn[j-1][2]=j+2*M-1;
					else 	  nn[j-1][2]=j-1;
				}
				else if(n==2)
				{	
					nn[j-1][0]=j-1;
					if(jj==L) nn[j-1][1]=j-2*(L-1)*M-1;
					else      nn[j-1][1]=j+2*M-1;
					if(ii==M) nn[j-1][2]=j-2*M+1;
					else 	  nn[j-1][2]=j+1;
				}
			}
		}
	}
	
	for(ii=0;ii<=N-1;ii++)
	{
		nnn[ii][0]=nn[nn[ii][0]-1][1];
		nnn[ii][1]=nn[nn[ii][0]-1][2];
		nnn[ii][2]=nn[nn[ii][1]-1][0];
		nnn[ii][3]=nn[nn[ii][1]-1][2];
		nnn[ii][4]=nn[nn[ii][2]-1][0];
		nnn[ii][5]=nn[nn[ii][2]-1][1];
	}
		
//	for(ii=0;ii<=N-1;ii++)
//	{
//		printf("%d %d %d\n",nn[ii][0],nn[ii][1],nn[ii][2]);
//	}
	
//	for(ii=0;ii<=N-1;ii++)
//	{
//		printf("%d %d %d %d %d %d\n",nnn[ii][0],nnn[ii][1],nnn[ii][2],nnn[ii][3],nnn[ii][4],nnn[ii][5]);
//	}
	
	kslop=sqrt(3)/3; 
	kb=4*sqrt(3)*Pi/(9*a);	
	for(kk=0;kk<M;kk++)
	{
		for(ii=0;ii<L;ii++)
		{
			k_x1[kk][ii]=(-M/2.+1.+kk)/M*b1[0]+(-L/2.+1.+ii)/L*b2[0];
			k_y1[kk][ii]=(-M/2.+1.+kk)/M*b1[1]+(-L/2.+1.+ii)/L*b2[1];
			if(k_y1[kk][ii]>-kslop*k_x1[kk][ii]+(L+0.5)/L*kb && k_x1[kk][ii]>0)
			{
				k_x1[kk][ii]=k_x1[kk][ii]-b1[0];
				k_y1[kk][ii]=k_y1[kk][ii]-b1[1];
			}	
			if(k_y1[kk][ii]>kslop*k_x1[kk][ii]+(L-0.5)/L*kb && k_x1[kk][ii]<=0)
			{
				k_x1[kk][ii]=k_x1[kk][ii]+b2[0];
				k_y1[kk][ii]=k_y1[kk][ii]+b2[1];
			}
			if(k_y1[kk][ii]<kslop*k_x1[kk][ii]-(L+0.5)/L*kb && k_x1[kk][ii]>0)
			{
				k_x1[kk][ii]=k_x1[kk][ii]-b2[0];
				k_y1[kk][ii]=k_y1[kk][ii]-b2[1];
			}	
			if(k_y1[kk][ii]<-kslop*k_x1[kk][ii]-(L-0.5)/L*kb && k_x1[kk][ii]<=0)
			{
				k_x1[kk][ii]=k_x1[kk][ii]+b1[0];
				k_y1[kk][ii]=k_y1[kk][ii]+b1[1];
			}
		}
	}
	ll=0;
	for(jj=0;jj<M;jj++)
	{
		for(ii=0;ii<L;ii++)
		{
			if(k_y1[ii][jj]>=kslop*k_x1[ii][jj]+(L-0.5)/L*kb && k_x1[ii][jj]<=0)
			{
				k_x1[ii][jj]=k_x1[ii][jj]+b2[0];
				k_y1[ii][jj]=k_y1[ii][jj]+b2[1];
			}
			if(k_y1[ii][jj]<=-kslop*k_x1[ii][jj]-(L-0.5)/L*kb && k_x1[ii][jj]<=0)
			{
				k_x1[ii][jj]=k_x1[ii][jj]+b1[0];
				k_y1[ii][jj]=k_y1[ii][jj]+b1[1];
			}
		k[ll][0]=k_x1[ii][jj];
		k[ll][1]=k_y1[ii][jj];
		ll=ll+1;
		}
	}
	
	
	for(ii=0;ii<N/2;ii++)
	{
		f_k[ii]=2*cos(sqrt(3)*a*k[ii][1])+4*cos(3*a/2*k[ii][0])*cos(sqrt(3)*a/2*k[ii][1]);
		omega2_1[ii]=2*gamma/m*(12+f_k[ii])-12*gamma/m*sqrt(3+f_k[ii]);
		omega2_2[ii]=2*gamma/m*(12+f_k[ii])+12*gamma/m*sqrt(3+f_k[ii]);
		omega_1[ii]=sqrt(omega2_1[ii]);
		omega_2[ii]=sqrt(omega2_2[ii]);
	}
	
	for(i=0;i<=N/2-1;i++)
	{
		for(j=i+1;j<N/2;j++)
		{
			if(omega_1[i]>=omega_1[j])		
			{
				koder[0][0]=k[i][0];omega_1oder=omega_1[i];
				koder[0][1]=k[i][1];
				k[i][0]=k[j][0];omega_1[i]=omega_1[j];
				k[i][1]=k[j][1];
				k[j][0]=koder[0][0];omega_1[j]=omega_1oder;
				k[j][1]=koder[0][1];
			}
			if(omega_2[i]>omega_2[j])		
			{
				omega_2oder=omega_2[i];
				omega_2[i]=omega_2[j];
				omega_2[j]=omega_2oder;
			}
		}
	}


	for(i=0;i<N/2;i++)
	{
		omega[i]=omega_1[i];
		Omega2[i]=omega[i]*omega[i];
		omega[i+N/2]=omega_2[i];
		Omega2[i+N/2]=omega_2[i]*omega_2[i];
	}
	
//	FILE *fp =NULL;
//	fp=fopen("omega.dat","w");
//	for(ii=0;ii<N;ii++)
//	{
//		fprintf(fp,"%d %f\n",ii+1,omega[ii]);
//	}
//	for(ii=0;ii<N;ii++)
//	{
//		printf("%d %f\n",ii+1,omega[ii]);
//	}

	for(i=0;i<N/2;i++)
	{
	num[i]=sin(k[i][0]*R[0][0]+k[i][1]*R[0][1])+sin(k[i][0]*R[1][0]+k[i][1]*R[1][1])+sin(k[i][0]*R[2][0]+k[i][1]*R[2][1]);
	den[i]=cos(k[i][0]*R[0][0]+k[i][1]*R[0][1])+cos(k[i][0]*R[1][0]+k[i][1]*R[1][1])+cos(k[i][0]*R[2][0]+k[i][1]*R[2][1]);
	phi[i]=atan(num[i]/den[i]);
	e[0][0][i]=-sqrt(2)/2*(cos(phi[i])+I*sin(phi[i]));
	e[0][1][i]=sqrt(2)/2*(cos(phi[i])+I*sin(phi[i]));
	e[1][0][i]=-sqrt(2)/2;
	e[1][1][i]=-sqrt(2)/2;
	}

	for(ii=0;ii<N/2;ii++)
	{
		for(jj=0;jj<N;jj++)
		{
		pp1[jj]=sqrt(m/(N/2))*(cos(k[ii][0]*r[jj][0]+k[ii][1]*r[jj][1])-I*sin(k[ii][0]*r[jj][0]+k[ii][1]*r[jj][1]))*e[jj%2][0][ii];
		pp2[jj]=sqrt(m/(N/2))*(cos(k[ii][0]*r[jj][0]+k[ii][1]*r[jj][1])-I*sin(k[ii][0]*r[jj][0]+k[ii][1]*r[jj][1]))*e[jj%2][1][ii];
		V[jj][ii]=pp1[jj];
		V[jj][N-1-ii]=pp2[jj];
		}
	}
	
//	FILE *fp =NULL;
//	fp=fopen("V.dat","w");
//	for(jj=0;jj<N;jj++)
//	{
//		for(ii=0;ii<N;ii++)
//		{
//			fprintf(fp,"%f+i%f ",creal(V[jj][ii]),cimag(V[jj][ii]));
////			fprintf(fp,"%f ",creal(V[jj][ii]));
////			fprintf(fp,"%f ",cimag(V[jj][ii]));
//		}
//		fprintf(fp,"\n");	
//	}
	
	for (rr=0;rr<ensemble;rr++)
	{
		P_total=0;
		srand(time(NULL));
		for (ii=0;ii<N_excited_modes;ii++)
		{
			rand_E[ii]=rand()/(RAND_MAX+1.0);	
			P_total=P_total+rand_E[ii];
			phi_r[ii]=rand()*2*Pi/(RAND_MAX+1.0);
	//		printf("%d  %f %f %f\n",ii+1,rand_E[ii],P_total,phi_r[ii]);
		}
		
		H_total=0;
		for (ii=0;ii<N;ii++)
		{
			if (ii>=min_excited_mode && ii<=max_excited_mode)
			{
				pE[ii-min_excited_mode]=rand_E[ii-min_excited_mode]/P_total;
				Eg[ii]=pE[ii-min_excited_mode]*Escaled;
				Qk0[ii]=sqrt(2*Eg[ii]/Omega2[ii])+I*0;
				Pk0[ii]=sqrt(2*Eg[ii])+I*0;
			}
			else
			{
				Eg[ii]=0;
				Qk0[ii]=0+I*0;
				Pk0[ii]=0+I*0;
			}
			H_total=H_total+Eg[ii];
	//		printf("%d  %f %f %f\n",ii+1,Eg[ii],H_total,Escaled);
		}
		
		for(jj=0;jj<N;jj++)
		{
			q0[jj]=0;
			p0[jj]=0;
			for (ii=0;ii<N;ii++)
			{
				q0[jj]=q0[jj]+creal(Qk0[ii])*creal(V[jj][ii]*(cos(phi_r[ii])+I*sin(phi_r[ii])));
				p0[jj]=p0[jj]-creal(Pk0[ii])*cimag(V[jj][ii]*(cos(phi_r[ii])+I*sin(phi_r[ii])));
			}
			qt0[jj]=q0[jj];
			pt0[jj]=p0[jj];		
		}
		
	
		
		for(kk=0;kk<20;kk++)
		{
			H_total=0;
			for(ii=0;ii<N;ii++)
			{
				Dz[0]=qt0[ii]-qt0[nn[ii][0]-1];Dz[1]=qt0[ii]-qt0[nn[ii][1]-1];Dz[2]=qt0[ii]-qt0[nn[ii][2]-1];
				Dz[3]=qt0[nnn[ii][0]-1]-qt0[nn[ii][0]-1];Dz[4]=qt0[nnn[ii][1]-1]-qt0[nn[ii][0]-1];
				Dz[5]=qt0[nnn[ii][2]-1]-qt0[nn[ii][1]-1];Dz[6]=qt0[nnn[ii][3]-1]-qt0[nn[ii][1]-1];
				Dz[7]=qt0[nnn[ii][4]-1]-qt0[nn[ii][2]-1];Dz[8]=qt0[nnn[ii][5]-1]-qt0[nn[ii][2]-1];
				U[ii]=potentials(Dz,parameter);
				T[ii]=0.5*pt0[ii]*pt0[ii];
				H_total=H_total+T[ii]+U[ii];
			}
			lambda=H_total/Escaled;
			for (ii=0;ii<N;ii++)
			{
				qt0[ii]=qt0[ii]/sqrt(sqrt(lambda));
				pt0[ii]=pt0[ii]/sqrt(sqrt(lambda));	
			}
		}
		epsilon_real=H_total/N;
		
	//	printf("%f\n",lambda);
	//	printf("%f\n",H_total);
	//	FILE *fp =NULL;
	//	fp=fopen("zvt.dat","w");
	//	for(jj=0;jj<N;jj++)
	//	{
	//		fprintf(fp,"%f\t",qt0[jj]);	
	//	}	
	//	fprintf(fp,"\n");	
	//	for(jj=0;jj<N;jj++)
	//	{
	//		fprintf(fp,"%f\t",pt0[jj]);	
	//	}
			
			for (ii=0;ii<N;ii++)
			{
				Dz[0]=qt0[ii]-qt0[nn[ii][0]-1];Dz[1]=qt0[ii]-qt0[nn[ii][1]-1];Dz[2]=qt0[ii]-qt0[nn[ii][2]-1];
				Dz[3]=qt0[nnn[ii][0]-1]-qt0[nn[ii][0]-1];Dz[4]=qt0[nnn[ii][1]-1]-qt0[nn[ii][0]-1];
				Dz[5]=qt0[nnn[ii][2]-1]-qt0[nn[ii][1]-1];Dz[6]=qt0[nnn[ii][3]-1]-qt0[nn[ii][1]-1];
				Dz[7]=qt0[nnn[ii][4]-1]-qt0[nn[ii][2]-1];Dz[8]=qt0[nnn[ii][5]-1]-qt0[nn[ii][2]-1];
				force[ii]=Forces(Dz,parameter);
				qt1[ii]=qt0[ii]+pt0[ii]+0.5*force[ii];
			}

			for (ii=0;ii<N;ii++)
			{
				Dz[0]=qt1[ii]-qt1[nn[ii][0]-1];Dz[1]=qt1[ii]-qt1[nn[ii][1]-1];Dz[2]=qt1[ii]-qt1[nn[ii][2]-1];
				Dz[3]=qt1[nnn[ii][0]-1]-qt1[nn[ii][0]-1];Dz[4]=qt1[nnn[ii][1]-1]-qt1[nn[ii][0]-1];
				Dz[5]=qt1[nnn[ii][2]-1]-qt1[nn[ii][1]-1];Dz[6]=qt1[nnn[ii][3]-1]-qt1[nn[ii][1]-1];
				Dz[7]=qt1[nnn[ii][4]-1]-qt1[nn[ii][2]-1];Dz[8]=qt1[nnn[ii][5]-1]-qt1[nn[ii][2]-1];
				force[ii]=Forces(Dz,parameter);
				qt2[ii]=2*qt1[ii]-qt0[ii]+force[ii];
				zvt[0][ii]=(qt2[ii]-qt0[ii])/2;
				zt[0][ii]=qt1[ii];	
			}
		
			for (jj=1;jj<checknum+1;jj++)
			{
				for(i=0;i<DN-1;i++)
				{
					for(ii=0;ii<N;ii++)	
					{
						qt0[ii]=qt1[ii];
						qt1[ii]=qt2[ii];
					}
					for (ii=0;ii<N;ii++)
					{
						Dz[0]=qt1[ii]-qt1[nn[ii][0]-1];Dz[1]=qt1[ii]-qt1[nn[ii][1]-1];Dz[2]=qt1[ii]-qt1[nn[ii][2]-1];
						Dz[3]=qt1[nnn[ii][0]-1]-qt1[nn[ii][0]-1];Dz[4]=qt1[nnn[ii][1]-1]-qt1[nn[ii][0]-1];
						Dz[5]=qt1[nnn[ii][2]-1]-qt1[nn[ii][1]-1];Dz[6]=qt1[nnn[ii][3]-1]-qt1[nn[ii][1]-1];
						Dz[7]=qt1[nnn[ii][4]-1]-qt1[nn[ii][2]-1];Dz[8]=qt1[nnn[ii][5]-1]-qt1[nn[ii][2]-1];
						force[ii]=Forces(Dz,parameter);
						qt2[ii]=2*qt1[ii]-qt0[ii]+force[ii];
					}			
				}			
				for(ii=0;ii<N;ii++)	
				{
					qt0[ii]=qt1[ii];
					qt1[ii]=qt2[ii];
				}
				for (ii=0;ii<N;ii++)
				{
					Dz[0]=qt1[ii]-qt1[nn[ii][0]-1];Dz[1]=qt1[ii]-qt1[nn[ii][1]-1];Dz[2]=qt1[ii]-qt1[nn[ii][2]-1];
					Dz[3]=qt1[nnn[ii][0]-1]-qt1[nn[ii][0]-1];Dz[4]=qt1[nnn[ii][1]-1]-qt1[nn[ii][0]-1];
					Dz[5]=qt1[nnn[ii][2]-1]-qt1[nn[ii][1]-1];Dz[6]=qt1[nnn[ii][3]-1]-qt1[nn[ii][1]-1];
					Dz[7]=qt1[nnn[ii][4]-1]-qt1[nn[ii][2]-1];Dz[8]=qt1[nnn[ii][5]-1]-qt1[nn[ii][2]-1];
					force[ii]=Forces(Dz,parameter);
					qt2[ii]=2*qt1[ii]-qt0[ii]+force[ii];
				}
				
				for (ii=0;ii<N;ii++)
				{
					zvt[jj][ii]=(qt2[ii]-qt0[ii])/2;
					zt[jj][ii]=qt1[ii];		
				}
			}
		
		for (ii=0;ii<checknum+1;ii++)
		{
			H_bar[ii]=0;
			for (jj=0;jj<N;jj++)
			{
				Pkt[ii][jj]=0;
				Qkt[ii][jj]=0;
				for (kk=0;kk<N;kk++)
				{
					Pkt[ii][jj]=Pkt[ii][jj]+zvt[ii][kk]*V[kk][jj];
					Qkt[ii][jj]=Qkt[ii][jj]+zt[ii][kk]*V[kk][jj];
				}
				HFt[ii][jj]=creal(0.5*Pkt[ii][jj]*conj(Pkt[ii][jj])+0.5*Omega2[jj]*Qkt[ii][jj]*conj(Qkt[ii][jj]));
				logHFt[ii][jj]=log10(HFt[ii][jj]);
				H_bar[ii]=H_bar[ii]+logHFt[ii][jj]/N;
			}
		}
		

		
		for (ii=1;ii<checknum+1;ii++)
		{
			SSigma1=0;
			for (jj=1;jj<N;jj++)
			{
				SSigma1=SSigma1+(logHFt[ii][jj]-H_bar[ii])*(logHFt[ii][jj]-H_bar[ii])/(N-1);
			}
			Sigma1[rr][ii]=sqrt(SSigma1);
		}
		
		for (ii=1;ii<checknum+1;ii++)
		{
			sumEktbar=0;
			kstart=floor(2*ii/3);
			Hkbar[ii]=0;
			for (jj=1;jj<N;jj++)
			{
				Ektbar[jj]=0;
				for (kk=kstart;kk<ii+1;kk++)
				{
					Ektbar[jj]=Ektbar[jj]+HFt[kk][jj]/(ii-kstart+1);		
				}
				sumEktbar=sumEktbar+Ektbar[jj];
					
//				logEktbar[ii][jj]=log10(Ektbar[jj]);
//				Hkbar[ii]=Hkbar[ii]+logEktbar[ii][jj]/N;								
			}
			
			eta=0;
			for (jj=1;jj<N;jj++)
			{
				Wk[jj]=Ektbar[jj]/sumEktbar;
				eta=eta-Wk[jj]*log(Wk[jj]);
			}
			Xi[rr][ii]=exp(eta)/N;		
		}
		
//		for (ii=1;ii<checknum+1;ii++)
//		{
//			SSigma2=0;
//			for (jj=1;jj<N;jj++)
//			{
//				SSigma2=SSigma2+(logEktbar[ii][jj]-Hkbar[ii])*(logEktbar[ii][jj]-Hkbar[ii])/(N-1);
//			}
//			Sigma2[rr][ii]=sqrt(SSigma2);
//		}
			
		printf("ensemble:%d\n",rr+1);
	}

	FILE *fp1 =NULL;
	fp1=fopen("HFt.dat","w");
	for (jj=0;jj<checknum+1;jj++)
	{
		for (ii=0;ii<N;ii++)
		{
			fprintf(fp1,"%f\t",HFt[jj][ii]);
		}	
		fprintf(fp1,"\n");
	}
	
	

	FILE *fp2 =NULL;
	fp2=fopen("Sigma1.dat","w");
	FILE *fp3 =NULL;
	fp3=fopen("Sigma2.dat","w");
	FILE *fp4 =NULL;
	fp4=fopen("Xi.dat","w");
	for (rr=0;rr<ensemble;rr++)
	{
		for (ii=0;ii<checknum+1;ii++)
		{
			fprintf(fp2,"%.15f\t",Sigma1[rr][ii]);
			fprintf(fp3,"%.15f\t",Sigma2[rr][ii]);
			fprintf(fp4,"%.15f\t",Xi[rr][ii]);
		}
		fprintf(fp2,"\n");
		fprintf(fp3,"\n");
		fprintf(fp4,"\n");
	}
	
	finish_run = clock();
	time_count=(double)(finish_run-start_run)/CLOCKS_PER_SEC;
	
	FILE *fp5 =NULL;
	fp5=fopen("parameter.dat","w");	
	fprintf(fp5,"%d\n%d\n%d\n%d\n%.6e\n%.6e\n%.6e\n%.4f\n",N,DN,checknum,iteration,t0,H_total,epsilon_real,time_count);	
	fclose(fp1);fclose(fp2);fclose(fp3);fclose(fp4);fclose(fp5);
	return 0;
}


