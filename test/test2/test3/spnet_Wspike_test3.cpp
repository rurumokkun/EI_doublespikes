// SPNET: Spiking neural network with axonal conduction delays and STDP
// Created by Eugene M. Izhikevich, May 17, 2004, San Diego, CA
// Saves spiking data each second in file spikes.dat
// To plot spikes, use MATLAB code: load spikes.dat;plot(spikes(:,1),spikes(:,2),'.');
#include <iostream>
#include <math.h>
#include <stdio.h>		
#include <sys/time.h>
#include <sys/resource.h>
using namespace std;
#include <stdlib.h> 
#include <time.h>
#include <cmath>

#define getrandom(max1) ((rand()%(int)((max1)))) // random integer between 0 and max-1     //getrandom(N)は0~N-1のいずれかを出力
#define NM 256
const   int		Ngr=2;                  	//ニューロングループの数 
const   int		SELF_SPIKE_TIME =1000;	                                //グループ別自己発火時間
const	int		INTERACT_TIME	=2000;					//グループ間相互作用時間
const	int		SIM_TIME=SELF_SPIKE_TIME+INTERACT_TIME;                //シミュレーション時間(シミュレーション前に決定)
const	int		NUM_EI_PATTERN = 10;	//				//EIの組み合わせの数  181031 10->1
const   int             num_trial=20;             //171101追記  1214に10->1へテスト用                     181031 20->1			 
const	int		N  = 1000;		// total number of neurons	
const   int     	interM=3;		// 1ニューロンあたりのnyu-ronnグループ間シナプス数		
const	int		M  = 100;		// the number of synapses per neuron 　
const	int		D  = 20;		// maximal axonal conduction delayst
		float	sm = 10.0;		// maximal synaptic strength		
int		post[Ngr][N][M];				// indeces of postsynaptic neurons
float	s[Ngr][N][M], sd[Ngr][N][M];		// matrix of synaptic weights and their derivatives
short	delays_length[Ngr][N][D];	// distribution of delays
short	delays[Ngr][N][D][M];		// arrangement of delays   
int		N_pre[Ngr][N], I_pre[Ngr][N][3*M], D_pre[Ngr][N][3*M];	// presynaptic information
int		interN_pre[Ngr][N], interI_pre[Ngr][N][10*interM], interD_pre[Ngr][N][10*interM];
float	*s_pre[Ngr][N][3*M], *sd_pre[Ngr][N][3*M];		// presynaptic weights
float	*inter_s_pre[Ngr][N][30*interM], *inter_sd_pre[Ngr][N][30*interM];		// グループ間シナプス前結合荷重&更新用変数
float	LTP[Ngr][N][1001+D], LTD[Ngr][N];	// STDP functions : LTP(Long Term Potentiation) LTD(Long Term Potentiation)
float	inter_LTP[Ngr][N][1001+D+10/*=30-10*/],inter_LTD[Ngr][N];		//	グループ間用のLTP,LTD
float	a[Ngr][N], d[Ngr][N];				// neuronal dynamics parameters
float	v[Ngr][N], u[Ngr][N];				// activity variables
int		N_firings[Ngr];				// the number of fired neurons 
const int N_firings_max=100*N/**10*/;	// upper limit on the number of fired neurons per sec		//発火数の上限: 変えていいのか？論文を確認すべき 9:1は必ず超える *10追加
int		firings[Ngr][N_firings_max][2]; // indeces and timings of spikes
int         interpost[Ngr][N][interM];		//
float       interpost_s[Ngr][N][interM];	//グループ間結合重み
float       interpost_sd[Ngr][N][interM];	//グループ間結合重みの更新用変数
float       interpost_delays[Ngr][N][interM];	//グループ間伝達遅れ
float       I[Ngr][N];				// 入力
int  		temp[Ngr]={0,0};			//
int		temp_gr;

float           full_LAP[Ngr][1000];//full_LAPはLAPデータの全時間分格納		//secごとにリセット
float		full_LAD[Ngr][1000];//抑制膜電流の平均
float           LAP[Ngr][1000];
float		LAD[Ngr][1000];

int		z[Ngr];				// N_firingsコピー用
int	EI_1,EI_2;
int	num;
int               Ne[Ngr][NUM_EI_PATTERN]={{800,800,800,800,800,800,800,800,800,800},     {800,800,800,800,800,800,800,800,800,800}};

int               Ni[Ngr][NUM_EI_PATTERN];





//////////////////////////////
//	grtoogr     //////////
//////////////////////////////
int grtoogr(int gr)
{
  int ogr;
   if(gr==0) {ogr=1;
    
   }
   else     { ogr=0;
    
   }
    return(ogr);
}



///////////////////////////////
///       初期化	///////
///////////////////////////////
void initialize()				
{	int i,j,k,jj,dd, exists,gr, r ; 		//gr:ニューロングループ指定用変数
        
         
	for(gr=0;gr<Ngr;gr++){
         for (i=0;i<Ne[gr][num];i++) a[gr][i]=0.02;// RS type
	for (i=Ne[gr][num];i<N;i++) a[gr][i]=0.1;  // FS type

	for (i=0;i<Ne[gr][num];i++) d[gr][i]=8.0;  // RS type
	for (i=Ne[gr][num];i<N;i++) d[gr][i]=2.0;  // FS type


			for(i=0;i<N;i++){              //181031				
                                LAP[gr][i]=0; //181031 shokika t[0 999] yori i[0 999] ga kanou
                                LAD[gr][i]=0;
                                full_LAP[gr][i]=0;
                                full_LAD[gr][i]=0;
                                
                         }
 

	srand((unsigned) time(NULL));
	for (i=0;i<N;i++) for (j=0;j<M;j++) 
	{
		do{
			exists = 0;		// avoid multiple synapses
			if (i<Ne[gr][num]) r = getrandom(N);
			else	  r = getrandom(Ne[gr][num]);// inh -> exc only
			if (r==i) exists=1;									// no self-synapses 
			for (k=0;k<j;k++) if (post[gr][i][k]==r) exists = 1;	// synapse already exists  
		}while (exists == 1);
		post[gr][i][j]=r;
	}
	for (i=0;i<Ne[gr][num];i++)	for (j=0;j<M;j++) s[gr][i][j]=6.0;  // initial exc. 興奮性シナプス荷重
	for (i=Ne[gr][num];i<N;i++)	for (j=0;j<M;j++) s[gr][i][j]=-(num+1); // 181015	yokusei　-1　to　-10	    抑制性シナプス荷重
  	for (i=0;i<N;i++)	for (j=0;j<M;j++) sd[gr][i][j]=0.0; // synaptic derivatives 
  	for (i=0;i<N;i++) 	     		//delays, delays_lengthの初期化:  興奮性,抑制性で場合分け
	{
		short ind=0;
		if (i<Ne[gr][num])
		{
			for (j=0;j<D;j++) 
			{	delays_length[gr][i][j]=M/D;	// uniform distribution of exc. synaptic delays
				for (k=0;k<delays_length[gr][i][j];k++)  // 0<=k<5
					delays[gr][i][j][k]=ind++;
			}
		}
		else
		{
			for (j=0;j<D;j++) delays_length[gr][i][j]=0;
			delays_length[gr][i][0]=M;			// all inhibitory delays are 1 ms
			for (k=0;k<delays_length[gr][i][0];k++)
					delays[gr][i][0][k]=ind++;
		}
	}
	
  	for (i=0;i<N;i++)
	{
		N_pre[gr][i]=0;
		for (j=0;j<Ne[gr][num];j++)
		for (k=0;k<M;k++)
		if (post[gr][j][k] == i)		// find all presynaptic neurons 
		{
			I_pre[gr][i][N_pre[gr][i]]=j;	// add this neuron to the list
			for (dd=0;dd<D;dd++)	// find the delay
				for (jj=0;jj<delays_length[gr][j][dd];jj++)
					if (post[gr][j][delays[gr][j][dd][jj]]==i) D_pre[gr][i][N_pre[gr][i]]=dd;	//0~19
			s_pre[gr][i][N_pre[gr][i]]=&s[gr][j][k];	// pointer to the synaptic weight	
			sd_pre[gr][i][N_pre[gr][i]++]=&sd[gr][j][k];// pointer to the derivative
		}
	}

	for (i=0;i<N;i++)	for (j=0;j<1+D;j++) LTP[gr][i][j]=0.0;   //伝達遅れが0~20それぞれに対して0で初期化
	for (i=0;i<N;i++)	LTD[gr][i]=0.0;
	for (i=0;i<N;i++)	v[gr][i]=-65.0;		// initial values for v
	for (i=0;i<N;i++)	u[gr][i]=0.2*v[gr][i];	// initial values for u

	N_firings[gr]=1;		// spike timings
	firings[gr][0][0]=-D-10;	// put a dummy spike at -D for simulation efficiency   181127 -D->-D+10
	firings[gr][0][1]=0;	// index of the dummy spike
																																																																																																																																																																																																																																																																																									
    }	//grについてのfor文終わり
  }   //initialize関数の終わり



////////////////////////////////
///	グループ初期化	////////
////////////////////////////////
void group_initialize()                 
{
  int gr,i,j,r,k,exists;

  
   for(gr=0;gr<Ngr;gr++){
    for(i=0;i<N;i++){
     for(j=0;j<interM;j++){				
       do{
         exists=0;
        if(i<Ne[gr][num])r=getrandom(N);				
        else    r=getrandom(Ne[gr][num]);
    	for(k=0;k<j;k++)if (interpost[gr][i][k]==r)exists=1;		
   
       }while(exists==1);
        interpost[gr][i][j]=r;					//シナプス後ニューロン
        if(i<Ne[gr][num])  interpost_s[gr][i][j]=6.0;				//重みの初期化(固定)
        else      interpost_s[gr][i][j]=-(1+num);
	interpost_sd[gr][i][j]=0;					//11/08追加
   	interpost_delays[gr][i][j]=getrandom(21)+10.0;			//伝達遅れの初期化: 10~30の間
      }//jに関するforループ
    }//iに関するforループ


   }//grに関するforループ


//////////////		グループ間の種々の変数の決定		/////////////////////////////////////////////////////////////
       for(gr=0;gr<Ngr;gr++)
	for(i=0;i<N;i++)
	{
		interN_pre[gr][i]=0;
		for(j=0;j<Ne[grtoogr(gr)][num];j++)
		for(k=0;k<interM;k++)
		if(interpost[grtoogr(gr)][j][k]==i)	//iはグループ(gr)
		{
			interI_pre[gr][i][interN_pre[gr][i]]=j;   //ニューロン //jはグループogr
		
			interD_pre[gr][i][interN_pre[gr][i]]=interpost_delays[grtoogr(gr)][j][k];	//10~30の間

			inter_s_pre[gr][i][interN_pre[gr][i]]=&interpost_s[grtoogr(gr)][j][k];
			inter_sd_pre[gr][i][interN_pre[gr][i]++]=&interpost_sd[grtoogr(gr)][j][k];
		}



	}


/////////////   2/1追加	/////////////////////////////////////////////////////////////////////////////////////////////
	for(gr=0;gr<Ngr;gr++)for (i=0;i<N;i++)	for (j=0;j<1+10+D;j++) inter_LTP[gr][i][j]=0.0;   //伝達遅れが10~30それぞれに対して0で初期化
	for(gr=0;gr<Ngr;gr++)for (i=0;i<N;i++)	inter_LTD[gr][i]=0.0;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        

}


int main()
{      
        char fname_spikes_1[50],fname_spikes_2[50];
        char fname_LAP_1[50],fname_LAP_2[50];
	char fname_interpost_s_1[50];
	char fname_interpost_s_2[50];	
	char fname_s_1[50];
	char fname_s_2[50];
       
        char str[30];
        char st[N]={'\0'};
        
	int		i, j, k, sec, t,gr,l,trial;
	float		p1,p2,p3,p4,q1,q2,q3,q4;
	time_t		t1,t2,timer;
        struct 		tm *timeptr;

        float		temp_I[Ngr][N];	
	float		v_sum[Ngr];
		


        
        int            firing_count[NUM_EI_PATTERN][Ngr][N];
        float          firing_rate_seq[NUM_EI_PATTERN][Ngr][SIM_TIME];  

	FILE	*fs_1,*fs_2,*fs_LAP_1,*fs_LAP_2;			      //fs_1はグループ１の発火情報, fs_2はグループ2の発火情報を出力する用
	FILE	*fp_interpost_s_1;			//グループ間結合重み
	FILE	*fp_interpost_s_2;			//グループ間結合重み
	FILE    *fp_s_1,*fp_s_2;					//グループ内結合重み
        FILE    *fp_firing_count;
	FILE    *fp_firing_rate;    //17'11/7



/////////////////////////////////////////////////////////////////////////////////////////// 	

	
	for(num=0;num<NUM_EI_PATTERN;num++){ 
	
	  Ni[0][num]=1000-Ne[0][num];
	  Ni[1][num]=1000-Ne[1][num];
	  
	}

//////////////////////////////////////////////////////////////////////////////////////////

	t1=time(NULL);		//開始時間

 for(trial=0;trial<num_trial;trial++){   //trial文
         timer=time(NULL);
         timeptr=localtime(&timer);
         strftime(st,NM,"data%Y%m%d%H%M",timeptr);

        system("mkdir data_s");
	system("mkdir data_interpost_s");
	system("mkdir data_LAP");
	system("mkdir data_spikes");
        system("mkdir data_firing_count");
	system("mkdir data`date '+%Y%m%d%H%M'`");
        sprintf(str,"mv data_* %s",st);
	 
 
     
     for(num=0;num<NUM_EI_PATTERN;num++){

        EI_1=Ne[0][num]*10+Ni[0][num]/10;
	EI_2=Ne[1][num]*10+Ni[1][num]/10;

        printf("EIバランスは,グループ1が%.1f:%.1f\tグループ2が%.1f:%.1fです\n",float(float(Ne[0][num])/float(100)),float(float(Ni[0][num])/float(100)),float(float(Ne[1][num])/float(100)),float(float(Ni[1][num])/float(100)));
	initialize();					// 関数呼び出し: assign connections, weights, etc. 
        group_initialize();				//　関数呼び出し: グループ初期化
 

        for(sec=0;sec<SIM_TIME;sec++){
  
          for(t=0;t<1000;t++){

           for(gr=0;gr<Ngr;gr++){
             for (i=0;i<N;i++) I[gr][i] = 0.0;	// reset the input 
			for (k=0;k<N/1000;k++)			// k? (k=0のみ?)
			  I[gr][getrandom(N)]=20.0;		// random thalamic :input 20の入力

            for(i=0;i<N;i++){
             if(v[gr][i]>=30){

               v[gr][i]=-65.0;

               u[gr][i]+=d[gr][i];

               LTP[gr][i][t+D]=0.1;

               LTD[gr][i]=0.12;
               inter_LTP[gr][i][t+D+10]=0.1;
           
               inter_LTD[gr][i]=0.12;


               for (j=0;j<N_pre[gr][i];j++) *sd_pre[gr][i][j]+=LTP[gr][I_pre[gr][i][j]][t+D-(D_pre[gr][i][j]+1)];// this spike was after pre-synaptic spikes    //sdの更新(pre->postの順に発火)

                         temp_gr=gr;
////////////////////////// updating inter_sd  /////////////////////////////////////////////////////////////////////////////////////////////////////
               if(sec>=SELF_SPIKE_TIME && gr==Ngr-1){
                 for(gr=0;gr<Ngr;gr++)
                  for(j=0;j<interN_pre[gr][i];j++)
                      *inter_sd_pre[gr][i][j]+=inter_LTP[grtoogr(gr)][interI_pre[gr][i][j]][t+D+10-interD_pre[gr][i][j]];
                 


                }
                          gr=temp_gr;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                      
                                firings[gr][N_firings[gr]  ][0]=t;			//[0]: 時間を格納
				firings[gr][N_firings[gr]++ ][1]=i;			//[1]: ニューロン番号を格納
						
				if (N_firings[gr] == N_firings_max) {cout << "Two many spikes at t=" << t << " (ignoring all)";N_firings[gr]=1;}
 



             }// if v[gr][i]>=30
           }// for(i)

                        z[gr]=N_firings[gr];
  			temp[gr]=z[gr];
			while (t-firings[gr][--z[gr]][0] <D) // jikan okure no jougen made sakanoboru
			{
				for (j=0; j< delays_length[gr][firings[gr][z[gr]][1]][t-firings[gr][z[gr]][0]]; j++)
				{
					i=post[gr][firings[gr][z[gr]][1]][delays[gr][firings[gr][z[gr]][1]][t-firings[gr][z[gr]][0]][j]]; 
					I[gr][i]+=s[gr][firings[gr][z[gr]][1]][delays[gr][firings[gr][z[gr]][1]][t-firings[gr][z[gr]][0]][j]];
					if (firings[gr][z[gr]][1] <Ne[gr][num]) // this spike is before postsynaptic spikes
						sd[gr][firings[gr][z[gr]][1]][delays[gr][firings[gr][z[gr]][1]][t-firings[gr][z[gr]][0]][j]]-=LTD[gr][i];   //sdの更新(post->preの順に発火)
					
					
				}
			}
                       temp_gr=gr;
                        
                      for(gr=0;gr<Ngr;gr++)   z[gr]=temp[gr];

                        gr=temp_gr;

/////////////////////////////////////////	interpost_sdの更新  	(021大きく修正 ) ////////////////////////////////////////// //0208コメントアウト
//	グループ内のsd更新の後、グループ間のsd更新		///	
		if(sec>=SELF_SPIKE_TIME && gr==1){
           
			  for(gr=0;gr<Ngr;gr++){
			  
				
  				while(t-firings[grtoogr(gr)][/*--*/z[grtoogr(gr)]][0] <D+10/* && z[grtoogr(gr)]!=0*/)  		                                                                                    
				{
					            
					if(t-firings[grtoogr(gr)][z[grtoogr(gr)]][0]>=10/* &&firings[grtoogr(gr)][z[grtoogr(gr)]][1]<Ne[grtoogr(gr)][num]*/){  //&&のあとz[]をfirings[][z[]][1]に変更 171214
					   for(j=0; j< interM; j++)
    					   {
	
                                               if(t-firings[grtoogr(gr)][z[grtoogr(gr)]][0]==interpost_delays[grtoogr(gr)][firings[grtoogr(gr)][z[grtoogr(gr)]][1]][j])
                                     		 {i=interpost[grtoogr(gr)][firings[grtoogr(gr)][z[grtoogr(gr)]][1]][j];	//グループgrの発火したニューロンfir..[1]のsyn_j先///shsei   
					            I[gr][i]+=interpost_s[grtoogr(gr)][firings[grtoogr(gr)][z[grtoogr(gr)]][1]][j];	//current must be  

						  if(firings[grtoogr(gr)][z[grtoogr(gr)]][1]<Ne[grtoogr(gr)][num])		//
					    	     interpost_sd[grtoogr(gr)][firings[grtoogr(gr)][z[grtoogr(gr)]][1]][j]-=inter_LTD[gr][i];     
                                                 }//if
				      		                                      /*}*/
				           }//for
                                        }//if
                                 z[grtoogr(gr)]--;
			       }//while



			  }//for
		

                          
			}//if		
                            
	       for(gr=0;gr<Ngr;gr++)   z[gr]=temp[gr];  ////
 				
		gr=temp_gr;	            
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                      
                  
                                v_sum[gr]=0;

				for(i=0;i<Ne[gr][num];i++){
                        	 v_sum[gr]+=v[gr][i];
				}
			
   				
				LAP[gr][t]=v_sum[gr]/float(Ne[gr][num]);
			        full_LAP[gr][t]=LAP[gr][t];

////////////////////////////////////// 抑制膜電流の平均 ////////////////////////////////////////////////////////////////////////


			        v_sum[gr]=0;

                                for(i=Ne[gr][num];i<N;i++){
                        	 v_sum[gr]+=v[gr][i];
				}
			
   				
				LAD[gr][t]=v_sum[gr]/float(Ni[gr][num]);
			        full_LAD[gr][t]=LAD[gr][t];
                
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////				


//////////////////////////    膜電位の更新(with オイラー法)       /////////////////////////////////////////////////////////////////////
///////////////////////////  time step is 0.5                     ////////////////////////////////////////////////


//                      if(gr==Ngr-1){  // I wo tukaunode 2tutomo no I ga sorowanto akan
//                       for(gr=0;gr<Ngr;gr++)
//		        for(i=0;i<N;i++){

////                         v[gr][i]+=0.5*((0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i]);
////                       
////                         v[gr][i]+=0.5*((0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i]);

//			k1=((0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i]);
//			k2=a[gr][i]*(0.2*v[gr][i]-u[gr][i]);
//			v[gr][i]+=k1*0.1;
//			u[gr][i]+=k2*0.1;
//			k1=((0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i]);
//			k2=a[gr][i]*(0.2*v[gr][i]-u[gr][i]);
//			v[gr][i]+=k1*0.1;
//			u[gr][i]+=k2*0.1;
//			k1=((0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i]);
//			k2=a[gr][i]*(0.2*v[gr][i]-u[gr][i]);
//			v[gr][i]+=k1*0.1;
//			u[gr][i]+=k2*0.1;
//			k1=((0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i]);
//			k2=a[gr][i]*(0.2*v[gr][i]-u[gr][i]);
//			v[gr][i]+=k1*0.1;
//			u[gr][i]+=k2*0.1;
//			k1=((0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i]);
//			k2=a[gr][i]*(0.2*v[gr][i]-u[gr][i]);
//			v[gr][i]+=k1*0.1;
//			u[gr][i]+=k2*0.1;
//			k1=((0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i]);
//			k2=a[gr][i]*(0.2*v[gr][i]-u[gr][i]);
//			v[gr][i]+=k1*0.1;
//			u[gr][i]+=k2*0.1;
//			k1=((0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i]);
//			k2=a[gr][i]*(0.2*v[gr][i]-u[gr][i]);
//			v[gr][i]+=k1*0.1;
//			u[gr][i]+=k2*0.1;
//			k1=((0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i]);
//			k2=a[gr][i]*(0.2*v[gr][i]-u[gr][i]);
//			v[gr][i]+=k1*0.1;
//			u[gr][i]+=k2*0.1;
//			k1=((0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i]);
//			k2=a[gr][i]*(0.2*v[gr][i]-u[gr][i]);
//			v[gr][i]+=k1*0.1;
//			u[gr][i]+=k2*0.1;
//			k1=((0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i]);
//			k2=a[gr][i]*(0.2*v[gr][i]-u[gr][i]);
//			v[gr][i]+=k1*0.1;
//			u[gr][i]+=k2*0.1;
//	



//			

////			 v[gr][i]+=0.25*((0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i]);
////                       
////                         v[gr][i]+=0.25*((0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i]);

////			 v[gr][i]+=0.25*((0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i]);
////                       
////                         v[gr][i]+=0.25*((0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i]);

////                          

//                     
//                         LTP[gr][i][t+D+1]=0.95*LTP[gr][i][t+D];
//                         LTD[gr][i]*=0.95;

//                             if(sec>=SELF_SPIKE_TIME){
//                               inter_LTP[gr][i][t+D+1+10]=0.95*inter_LTP[gr][i][t+D+10];
//                               inter_LTD[gr][i]*=0.95;
//                             }



//                        }
//                        

//                       }

////                       gr=temp_gr;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////    膜電位の更新(with ルンゲクッタ法)       /////////////////////////////////////////////////////////////////////
///////////////////////  time step is 0.5                     ////////////////////////////////////////////////


                      if(gr==Ngr-1){  // I wo tukaunode 2tutomo no I ga sorowanto akan
                       for(gr=0;gr<Ngr;gr++)
		        for(i=0;i<N;i++){
			p1=(0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i];
			q1=a[gr][i]*(0.2*v[gr][i]-u[gr][i]);
			p2=(0.04*(v[gr][i]+1/2*p1*0.5)+5)*(v[gr][i]+1/2*p1*0.5)+140-(u[gr][i]+1/2*q1*0.5)+I[gr][i];
			q2=a[gr][i]*(0.2*(v[gr][i]+1/2*p1*0.5)-(u[gr][i]+1/2*q1*0.5));
			p3=(0.04*(v[gr][i]+1/2*p2*0.5)+5)*(v[gr][i]+1/2*p2*0.5)+140-(u[gr][i]+1/2*q2*0.5)+I[gr][i];
			q3=a[gr][i]*(0.2*(v[gr][i]+1/2*p2*0.5)-(u[gr][i]+1/2*q2*0.5));
			p4=(0.04*(v[gr][i]+p3*0.5)+5)*(v[gr][i]+p3*0.5)+140-(u[gr][i]+q3*0.5)+I[gr][i];
			q4=a[gr][i]*(0.2*(v[gr][i]+p3*0.5)-(u[gr][i]+q3*0.5));

			v[gr][i]+=0.5/6*(p1+2*p2+2*p3+p4);
			u[gr][i]+=0.5/6*(q1+2*q2+2*q3+q4);

			p1=(0.04*v[gr][i]+5)*v[gr][i]+140-u[gr][i]+I[gr][i];
			q1=a[gr][i]*(0.2*v[gr][i]-u[gr][i]);
			p2=(0.04*(v[gr][i]+1/2*p1*0.5)+5)*(v[gr][i]+1/2*p1*0.5)+140-(u[gr][i]+1/2*q1*0.5)+I[gr][i];
			q2=a[gr][i]*(0.2*(v[gr][i]+1/2*p1*0.5)-(u[gr][i]+1/2*q1*0.5));
			p3=(0.04*(v[gr][i]+1/2*p2*0.5)+5)*(v[gr][i]+1/2*p2*0.5)+140-(u[gr][i]+1/2*q2*0.5)+I[gr][i];
			q3=a[gr][i]*(0.2*(v[gr][i]+1/2*p2*0.5)-(u[gr][i]+1/2*q2*0.5));
			p4=(0.04*(v[gr][i]+p3*0.5)+5)*(v[gr][i]+p3*0.5)+140-(u[gr][i]+q3*0.5)+I[gr][i];
			q4=a[gr][i]*(0.2*(v[gr][i]+p3*0.5)-(u[gr][i]+q3*0.5));

			v[gr][i]+=0.5/6*(p1+2*p2+2*p3+p4);
			u[gr][i]+=0.5/6*(q1+2*q2+2*q3+q4);
			
 
                         LTP[gr][i][t+D+1]=0.95*LTP[gr][i][t+D];
                         LTD[gr][i]*=0.95;

                             if(sec>=SELF_SPIKE_TIME){
                               inter_LTP[gr][i][t+D+1+10]=0.95*inter_LTP[gr][i][t+D+10];
                               inter_LTD[gr][i]*=0.95;
                             }



                        }
                        

                       }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////



           }// for gr




          }// for t 
         
        
//////////////////////////        prepare for the next step    /////////////////////////////////////////////
          for(gr=0;gr<Ngr;gr++){

 
           cout << "sec=" << sec << ", firing rate=" << float(N_firings[gr])/N << "\tgroup1=   " << Ne[0][num]/100 << ":" << Ni[0][num]/100 << "\tgroup2=   " << Ne[1][num]/100 << ":" << Ni[1][num]/100 << "\t" << "trial#" << trial+1 << "\n";
           firing_rate_seq[num][gr][sec]=float(N_firings[gr])/N;

		   if(gr==0){
		    sprintf(fname_spikes_1,"spikes_1_%dsec_%d_%d.txt",SIM_TIME,EI_1,EI_2);

		    fs_1=fopen(fname_spikes_1,"w");

		    for(i=1;i<N_firings[gr];i++)
		      if(firings[gr][i][0]>=0)
		           fprintf(fs_1,"%d  %d\n", firings[gr][i][0],firings[gr][i][1]);
		    fclose(fs_1);
		   }
         

                   else{
		    sprintf(fname_spikes_2,"spikes_2_%dsec_%d_%d.txt",SIM_TIME,EI_1,EI_2);

		    fs_2=fopen(fname_spikes_2,"w");

		    for(i=1;i<N_firings[gr];i++)
		      if(firings[gr][i][0]>=0)
		           fprintf(fs_2,"%d  %d\n", firings[gr][i][0],firings[gr][i][1]);
		    fclose(fs_2);
		   }




                     for(i=0;i<N;i++)          // prepare for the next sec
                        for(j=0;j<D+1;j++)
			   LTP[gr][i][j]=LTP[gr][i][1000+j];


                     z[gr]=N_firings[gr]-1;
         
                     for(i=0;i<N;i++)
                        for(j=0;j<D+10+1;j++)
                           inter_LTP[gr][i][j]=inter_LTP[gr][i][1000+j];
                        
                     while(1000-firings[gr][z[gr]][0]<D+10/* 181109 D->D+10  181126 ->D*/) z[gr]--;
                     for(i=1;i<N_firings[gr]-z[gr];i++){                  
                          firings[gr][i][0]=firings[gr][z[gr]+i][0]-1000;
			  firings[gr][i][1]=firings[gr][z[gr]+i][1];
                     
                      } 


                       N_firings[gr]=N_firings[gr]-z[gr]; // (20msec sakanoboru baai)

                     if(sec<SIM_TIME-100){
                       for(i=0;i<Ne[gr][num];i++)
                        for(j=0;j<M;j++){
                         s[gr][i][j]+=0.01+sd[gr][i][j];
                         sd[gr][i][j]*=0.9;
                         if (s[gr][i][j]>sm) s[gr][i][j]=sm;
                         if (s[gr][i][j]<0)  s[gr][i][j]=0.0;


                        }

  
                     }


//////////////////////////// updating connection weights between two groups ///////////////////////////

                     if(sec>=SELF_SPIKE_TIME && sec<SIM_TIME-100){
                       for(i=0;i<Ne[gr][num];i++)
                         for(j=0;j<interM;j++){
                          interpost_s[gr][i][j]+=0.01+interpost_sd[gr][i][j];
                          interpost_sd[gr][i][j]*=0.9;
                          if(interpost_s[gr][i][j]>sm) interpost_s[gr][i][j]=sm;
                          if(interpost_s[gr][i][j]<0)  interpost_s[gr][i][j]=0.0;
                         
                         }


                     
                     }

           }
///////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////
///////             LAP&LAD                    ///////////
//////////////////////////////////////////////////////////
         sprintf(fname_LAP_1,"LAP_1_%d_%d_%d_%dsec.txt",EI_1,EI_2,-(num+1),SIM_TIME);
	 sprintf(fname_LAP_2,"LAP_2_%d_%d_%d_%dsec.txt",EI_1,EI_2,-(num+1),SIM_TIME);
	 
         fs_LAP_1=fopen(fname_LAP_1,"a");
	 fs_LAP_2=fopen(fname_LAP_2,"a");


              for(t=0;t<1000;t++){

                 fprintf(fs_LAP_1,"%f\t%f\n",full_LAP[0][t],full_LAD[0][t]);
                 fprintf(fs_LAP_2,"%f\t%f\n",full_LAP[1][t],full_LAD[1][t]);


              }
              fclose(fs_LAP_1);
              fclose(fs_LAP_2);


///////////////////////////////////////////////////////////

         if(sec==0||sec==1000||sec==1999||sec==2999){
       
           sprintf(fname_interpost_s_1,"interpost_s_%d_%dsec_%d_%d_%d.txt",1,sec,EI_1,EI_2,-(num+1));
	   sprintf(fname_interpost_s_2,"interpost_s_%d_%dsec_%d_%d_%d.txt",2,sec,EI_1,EI_2,-(num+1));
           sprintf(fname_s_1,"s_%d_%dsec_%d_%d_%d.txt",1,sec,EI_1,EI_2,-(num+1));
	   sprintf(fname_s_2,"s_%d_%dsec_%d_%d_%d.txt",2,sec,EI_1,EI_2,-(num+1));

           fp_interpost_s_1=fopen(fname_interpost_s_1,"w");
	   fp_interpost_s_2=fopen(fname_interpost_s_2,"w");
	   fp_s_1=fopen(fname_s_1,"w");
           fp_s_2=fopen(fname_s_2,"w");

             for(i=0;i<N;i++){
               for(j=0;j<interM;j++){
                   fprintf(fp_interpost_s_1,"%f\t%d\t%d\n",interpost_s[0][i][j],i,interpost[0][i][j]);

                }
              }
 
           
              for(i=0;i<N;i++){
               for(j=0;j<interM;j++){
                   fprintf(fp_interpost_s_2,"%f\t%d\t%d\n",interpost_s[Ngr-1][i][j],i,interpost[Ngr-1][i][j]);

                }
              }


             for(i=0;i<N;i++){
               for(j=0;j<M;j++){
                 fprintf(fp_s_1,"%f\t%d\t%d\n",s[0][i][j],i,post[0][i][j]);

               }

              }


              for(i=0;i<N;i++){
               for(j=0;j<M;j++){
                 fprintf(fp_s_2,"%f\t%d\t%d\n",s[Ngr-1][i][j],i,post[Ngr-1][i][j]);

               }

              }


              fclose(fp_interpost_s_1);
	      fclose(fp_interpost_s_2);
              fclose(fp_s_1);
              fclose(fp_s_2);




           }// for if

///////////////////////////////////////////////////////////////////


        }// for sec



        if(num<NUM_EI_PATTERN-1){
          cout<<"run the following simulation:"<<endl;

           for(i=0;i<Ngr;i++){
            printf("%d %d %d\n",Ne[i][num+1],Ni[i][num+1],-(num+2));
           }

        }

        else  cout<<"simulation has finished!!"<<endl;



    }// for num

     fp_firing_count=fopen("firing_count.txt","w");
     for(gr=0;gr<Ngr;gr++)
       for(i=0;i<N;i++)
         for(num=0;num<NUM_EI_PATTERN;num++)

           if(num<NUM_EI_PATTERN-1)fprintf(fp_firing_count,"%d\t",firing_count[num][gr][i]);
           else                    fprintf(fp_firing_count,"%d\n",firing_count[num][gr][i]);
           if(i==N-1)              fprintf(fp_firing_count,"\n\n");

           fclose(fp_firing_count);


     fp_firing_rate=fopen("firing_rate_seq.txt","w");
      for(gr=0;gr<Ngr;gr++)
       for(sec=0;sec<SIM_TIME;sec++)
        for(num=0;num<NUM_EI_PATTERN;num++)

            if(num<NUM_EI_PATTERN-1) 	fprintf(fp_firing_rate,"%f\t",firing_rate_seq[num][gr][sec]);
            else 			       fprintf(fp_firing_rate,"%f\n",firing_rate_seq[num][gr][sec]);



            fclose(fp_firing_rate);

        system("mv spikes* data_spikes"); 
	system("mv LAP* data_LAP");
	system("mv s_* data_s");
	system("mv interpost_s* data_interpost_s");
        system("mv firing_count.txt data_firing_count");
        system("mv firing_rate_seq.txt data_firing_count");
	system(str); 



 }// for trial
   t2=time(NULL);
   printf("time=%d[s]\n",(int)(t2-t1));


}// int main
