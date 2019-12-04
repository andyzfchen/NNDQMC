#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <string>
#include <Eigen/Eigen>
#include <random>
#include <unsupported/Eigen/MatrixFunctions>
#include <chrono>
#include <fstream>

using namespace std;
#define Dmatrix Eigen::MatrixXd
#define map Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >

int getparams(int argc, char *argv[], int &N, int &L, int &warmup, int &count, int &avgsteps, double &g, double &w, double &beta, int &lagmax){
	if(argc!=10){
		std::cout<<"ERROR!! input parameters in this order: lattice size, imaginary time slices, warm up steps, gcount, avgsteps, coupling strength, frequency, inverse temp, lagmax"<<std::endl;
		return 0;
	}
	std::istringstream Ns(string(argv[1])+" ");
	if(!(Ns>>N)){
		std::cout<<"ERROR with input N"<<std::endl;
		return 0;
	}
	std::istringstream Ls(string(argv[2])+" ");
	if(!(Ls>>L)){
		std::cout<<"ERROR with input L"<<std::endl;
		return 0;
	}
	std::istringstream wus(string(argv[3])+" ");
	if(!(wus>>warmup)){
		std::cout<<"ERROR with input warmup"<<std::endl;
		return 0;
	}
	std::istringstream counts(string(argv[4])+" ");
	if(!(counts>>count)){
		std::cout<<"ERROR with input count"<<std::endl;
		return 0;
	}
	std::istringstream avgs(string(argv[5])+" ");
	if(!(avgs>>avgsteps)){
		std::cout<<"ERROR with input avgsteps"<<std::endl;
		return 0;
	}
	std::istringstream gs(string(argv[6])+" ");
	if(!(gs>>g)){
		std::cout<<"ERROR with input g"<<std::endl;
		return 0;
	}
	std::istringstream ws(string(argv[7])+" ");
	if(!(ws>>w)){
		std::cout<<"ERROR with input w"<<std::endl;
		return 0;
	}
	std::istringstream bs(string(argv[8])+" ");
	if(!(bs>>beta)){
		std::cout<<"ERROR with input b"<<std::endl;
		return 0;
	}
	std::istringstream ls(string(argv[9])+" ");
	if(!(ls>>lagmax)){
		std::cout<<"ERROR with input lagmax"<<std::endl;
		return 0;
	}
	return 1;
}




// This function does the numerically stable calculation of the Green's function described in Appendix B of the writeup
double buildG(Dmatrix &G, double *B, int L, int N, int mmax, int l, int &lcount, double &avgerr, unsigned int &avgerrcounter){

	//break apart the multiplication of L B matricies into L/m multiplications of m matricies
	int lvls=L/mmax;
	int rem=L-mmax*lvls;
	int N4=N*N*N*N;
	int N2=N*N;


	//counters
	int lvl,m,i;


	//define A, D and R.
	Dmatrix A=Dmatrix::Identity(N2,N2);
	Eigen::VectorXd D=Eigen::VectorXd::Ones(N2);
	Dmatrix R=Dmatrix::Identity(N2,N2);
	Eigen::HouseholderQR<Dmatrix> qr;
	Dmatrix Rp(N2,N2);



	//at each level do the following: multiply m Bs to the current A
	//                                multiply the previous D to A
	//                                QR decompose A
	//                                set A to the orthogonal part, D to the diagonal of R, and R to rescaled R*old R
	int currm=l-1;
	for(lvl=0;lvl<lvls;lvl++){
		for(m=0;m<mmax;m++){
			A=map(B+(currm+L)%L*N4,N2,N2)*A;
			currm--;
		}
		A=A*D.asDiagonal();
		qr = A.householderQr();
		A=qr.householderQ()*Dmatrix::Identity(N2,N2);
		Rp=qr.matrixQR().triangularView<Eigen::Upper>();
		D=Rp.diagonal();
		for(i=0;i<N2;i++){
			Rp.row(i)/=D(i);
		}
		R=Rp*R;
	}


	//do the things described above on the remaining B matricies
	for(m=0;m<rem;m++){
		A=map(B+(currm+L)%L*N4,N2,N2)*A;
		currm--;
	}
	A=A*D.asDiagonal();
	qr = A.householderQr();
	A=qr.householderQ()*Dmatrix::Identity(N2,N2);
	Rp=qr.matrixQR().triangularView<Eigen::Upper>();
	D=Rp.diagonal();
	for(i=0;i<N2;i++){
		Rp.row(i)/=D(i);
	}
	R=Rp*R;


	//F=U^{-1}*R^{-1}+D
	//A= old A * U from F
	//R= R from F* old R
	//element wise inverse of D
	Dmatrix F=A.transpose()*R.triangularView<Eigen::Upper>().solve(Dmatrix::Identity(N2,N2))+(D.asDiagonal()*Dmatrix::Identity(N2,N2));
	qr=F.householderQr();
	A=A*qr.householderQ()*Dmatrix::Identity(N2,N2);
	Rp=qr.matrixQR().triangularView<Eigen::Upper>();
	D=Rp.diagonal();
	for(i=0;i<N2;i++){
		Rp.row(i)/=D(i);
		D(i)=1./D(i);
	}
	R=Rp*R;


	//put final inverse calculation into G
	F=R.triangularView<Eigen::Upper>().solve(Dmatrix::Identity(N2,N2))*D.asDiagonal()*A.transpose();
	double err=(G-F).norm()/N2/N2;
	G=F;
	lcount=0;
	avgerr=(avgerr*avgerrcounter+err)/(avgerrcounter+1.);
	avgerrcounter++;
	return err;
}





// This function performs the local update sweeps
void LocalUpdates(int &total, int &accept, double &errG, int &lcount, int m, int &L, int &N, double &g, double &dt, double &w, std::uniform_real_distribution<double> &unif, mt19937_64 &rng, double *X, Dmatrix &G, double *B, int &lcountmax, int &lcountb12, double &avgerr, unsigned int &avgerrcounter, Eigen::MatrixXd &WL1, Eigen::MatrixXd &WL2, Eigen::VectorXd &WL3, double &WLz, Eigen::VectorXd &bL1, Eigen::VectorXd &bL2, double &bL3, Eigen::VectorXd &xmean, Eigen::VectorXd &xstd, double &zmean, double &zstd, double &deaccum){

	//counters
	int l,yy,xx,loc;
	int N2=N*N;
	int N4=N*N*N*N;

	//useful values
	double dx,delta,prob;
	Eigen::VectorXd Input(28);
	Eigen::VectorXd First(15);
	Eigen::VectorXd Second(5);
	double out;
	double avgx=-3.9983295349324486;
	double z;

	//iterate through time slices
	for(l=0;l<L;l++){
		lcount++;
		loc=-1;
		//Iterate through the spatial sites
		for(yy=0;yy<N;yy++){
			for(xx=0;xx<N;xx++){

				//suggest move
				loc++;
				total++;
				dx=2.*(unif(rng)-0.5);
				delta=exp(-dt*g*dx)-1;
				//build Input vactor
				Input(0)=dx;
				Input(1)=X[l*N2+yy*N+xx];
				Input(2)=X[l*N2+yy*N+(xx+1)%N];
				Input(3)=X[l*N2+yy*N+(xx-1+N)%N];
				Input(4)=X[l*N2+((yy+1)%N)*N+xx];
				Input(5)=X[l*N2+((yy-1+N)%N)*N+xx];
				Input(6)=X[l*N2+((yy+1)%N)*N+(xx+1)%N];
				Input(7)=X[l*N2+((yy+1)%N)*N+(xx-1+N)%N];
				Input(8)=X[l*N2+((yy-1+N)%N)*N+(xx+1)%N];
				Input(9)=X[l*N2+((yy-1+N)%N)*N+(xx-1+N)%N];
				Input(10)=X[((l+1)%L)*N2+yy*N+xx];
				Input(11)=X[((l+1)%L)*N2+yy*N+(xx+1)%N];
				Input(12)=X[((l+1)%L)*N2+yy*N+(xx-1+N)%N];
				Input(13)=X[((l+1)%L)*N2+((yy+1)%N)*N+xx];
				Input(14)=X[((l+1)%L)*N2+((yy-1+N)%N)*N+xx];
				Input(15)=X[((l+1)%L)*N2+((yy+1)%N)*N+(xx+1)%N];
				Input(16)=X[((l+1)%L)*N2+((yy+1)%N)*N+(xx-1+N)%N];
				Input(17)=X[((l+1)%L)*N2+((yy-1+N)%N)*N+(xx+1)%N];
				Input(18)=X[((l+1)%L)*N2+((yy-1+N)%N)*N+(xx-1+N)%N];
				Input(19)=X[((l-1+L)%L)*N2+yy*N+xx];
				Input(20)=X[((l-1+L)%L)*N2+yy*N+(xx+1)%N];
				Input(21)=X[((l-1+L)%L)*N2+yy*N+(xx-1+N)%N];
				Input(22)=X[((l-1+L)%L)*N2+((yy+1)%N)*N+xx];
				Input(23)=X[((l-1+L)%L)*N2+((yy-1+N)%N)*N+xx];
				Input(24)=X[((l-1+L)%L)*N2+((yy+1)%N)*N+(xx+1)%N];
				Input(25)=X[((l-1+L)%L)*N2+((yy+1)%N)*N+(xx-1+N)%N];
				Input(26)=X[((l-1+L)%L)*N2+((yy-1+N)%N)*N+(xx+1)%N];
				Input(27)=X[((l-1+L)%L)*N2+((yy-1+N)%N)*N+(xx-1+N)%N];
				double c=X[l*N2+yy*N+xx]-avgx;
				double cp=c+dx;
				z=0.25*cp*cp*cp*cp-0.5*avgx*avgx*cp*cp-0.25*c*c*c*c+0.5*avgx*avgx*c*c;

				//scale inputs
				Input=((Input.array()-xmean.array())/xstd.array()).matrix();
				z=(z-zmean)/zstd;
				

				//calculate output
				First=(Input.adjoint()*WL1+bL1.adjoint()).adjoint();
				First=Eigen::tanh(First.array());
				Second=(First.adjoint()*WL2+bL2.adjoint()).adjoint();
				Second=Eigen::tanh(Second.array());
				out=Second.adjoint()*WL3+bL3+WLz;
				prob=min(1.,exp(-L*dt*out));

				//change B if accepted
				if(unif(rng)<prob){
					accept++;
					//Bp=(I+Delta)*B
					double *Bp=B+l*N4+loc*N2;
					for(int k=0;k<N2;k++){
						Bp[k]*=(1+delta);
					}
					X[l*N2+loc]+=dx;
					deaccum+=out;
				}
			}
		}
	}
}


//This function performs the global update sweeps
void GlobalUpdates(int &gtotal, int &gaccept, int &count, int lcount, int &m, int &L, int &N, double &g, double &dt, double &w, std::uniform_real_distribution<double> &unif, mt19937_64 &rng, double *X, Dmatrix &G, double *B, double *Wconv,	Eigen::MatrixXd &WG1, Eigen::VectorXd &WG2, double &WGg, double &WGz, Eigen::VectorXd &theta, Eigen::VectorXd &bconv, Eigen::VectorXd &bG1, double &bG2, Eigen::MatrixXd &xGmean, Eigen::MatrixXd &xGstd, double &dxGmean, double &dxGstd, double &zGmean, double &zGstd, double &gGmean, double &gGstd, double &deaccum){

	int N2=N*N;
	int N4=N*N*N*N;

	//counters
	int i,l,k,f;

	//values
	int sitex,sitey;
	double delta, *Bp, prob,c,cp;
	Eigen::VectorXd flat(120);
	Eigen::VectorXd flat2(20);
	double out,dx,z,gg,avgx=-3.9996732015947853;
	Eigen::MatrixXd Input(L,3);
	Eigen::MatrixXd Filter(8,3);

	//perform updates at count numer of sites
	for(i=0;i<count;i++){
		gtotal++;

		//the site we are doing the change on
		sitex=(int)floor(N*unif(rng));
		sitey=(int)floor(N*unif(rng));
		
		c=cp=g=z=0;
		//calculate g and z
		for(l=0;l<L;l++){
			c=X[l*N2+sitey*N+sitex]-avgx;
			cp=X[l*N2+sitey*N+sitex]-avgx+dx;
			gg+=0.25*cp*cp*cp*cp-0.5*avgx*avgx*cp*cp-0.25*c*c*c*c+0.5*avgx*avgx*c*c;
			z+=(cp+avgx)*(cp+avgx)-X[l*N2+sitey*N+sitex]*X[l*N2+sitey*N+sitex];
		}
		gg/=L;
		
		//calculate Input and scale the inputs
		dx=(unif(rng)-0.5)*2*g/w/w;
		for(l=0;l<L;l++){
			Input(l,0)=X[l*N2+sitey*N+sitex];
			Input(l,1)=0.25*(X[l*N2+sitey*N+(sitex+1)%N]+X[l*N2+sitey*N+(sitex-1+N)%N]+X[l*N2+((sitey+1)%N)*N+sitex]+X[l*N2+((sitey-1+N)%N)*N+sitex]);
			Input(l,2)=0.25*(X[l*N2+((sitey+1)%N)*N+(sitex+1)%N]+X[l*N2+((sitey+1)%N)*N+(sitex-1+N)%N]+X[l*N2+((sitey-1+N)%N)*N+(sitex+1)%N]+X[l*N2+((sitey-1+N)%N)*N+(sitex-1+N)%N]);
		}
		Input=((Input.array()-xGmean.array())/xGstd.array()).matrix();
		dx=(dx-dxGmean)/dxGstd;
		gg=(gg-gGmean)/gGstd;
		z=(z-zGmean)/zGstd;



		for(f=0;f<10;f++){
			Filter=map(Wconv,8,3);
			for(k=0;k<12;k++){
				flat(k*10+f)=(Filter.array()*Input.block(k*3,0,8,3).array()).sum()+bconv(f)+dx*theta(f);
			}
		}
		flat=Eigen::tanh(flat.array());
		flat2=(flat.adjoint()*WG1+bG1.adjoint()).adjoint();
		flat2=Eigen::tanh(flat2.array());
		out=flat2.adjoint()*WG2+bG2+WGz*z+WGg*gg;







		//either accept or reject
		prob=min(1.,exp(-L*dt*out));
		if(unif(rng)<prob){
			gaccept++;
			deaccum+=out;	
			for(l=0;l<L;l++){
				X[l*N2+sitey*N+sitex]+=dx;
			}
			delta=exp(-dt*g*dx)-1;
			for(l=0;l<L;l++){
				Bp=B+(l*N4+(sitey*N+sitex)*N2);
				for(k=0;k<N2;k++){
					Bp[k]*=(1+delta);
				}
			}
		}
	}
}


double calcS(double *X, double w, double dt, int N2, int L){
	double S=0;
	for(int l=0;l<L;l++){
		for(int n=0;n<N2;n++){
			S+=0.5*w*w*X[l*N2+n]*X[l*N2+n]+0.5/dt/dt*(X[l*N2+n]-X[((l+1)%L)*N2+n])*(X[l*N2+n]-X[((l+1)%L)*N2+n]);
		}
	}
	return S;
}

// This function performs the local update sweeps
void LocalUpdates(int &total, int &accept, double &errG, int &lcount, int m, int &L, int &N, double &g, double &dt, double &w, std::uniform_real_distribution<double> &unif, mt19937_64 &rng, double *X, Dmatrix &G, double *B, int &lcountmax, int &lcountb12, double &avgerr, unsigned int &avgerrcounter){

  //counters
  int l,yy,xx,loc;
  int N2=N*N;
  int N4=N*N*N*N;

  //useful values
  double dx, x, xm, xp,dS,delta,prob,opip;
  Dmatrix F(N2,N2);
  Dmatrix R(N2,N2);
  Dmatrix uvT(N2,N2);
  Eigen::VectorXd D(N2);
  Eigen::HouseholderQR<Dmatrix> qr;

  //iterate through time slices
  for(l=0;l<L;l++){
    lcount++;
    loc=-1;
    //Iterate through the spatial sites
    for(yy=0;yy<N;yy++){
      for(xx=0;xx<N;xx++){

        //suggest move
        loc++;
        total++;
        dx=2.*(unif(rng)-0.5);
        x=X[l*N2+loc];
        xp=X[(l+1)%L*N2+loc];
        xm=X[(l-1+L)%L*N2+loc];
        dS=1./2*w*w*dx*(2.*x+dx)+dx/dt/dt*(dx-xp-xm+2.*x);
        delta=exp(-dt*g*dx)-1;
        opip=1+delta*(1-G(loc,loc));
        prob=min(1.,exp(-dS*dt)*opip*opip);

        //change B and G if accepted
        if(unif(rng)<prob){
          accept++;
          //Bp=(I+Delta)*B
          double *Bp=B+l*N4+loc*N2;
          for(int k=0;k<N2;k++){
            Bp[k]*=(1+delta);
					}
          X[l*N2+loc]+=dx;
          //Gp=G(1-1/(1+vTu)uvT)
          uvT = Dmatrix::Zero(N2,N2);
          uvT.row(loc)=-delta*G.row(loc);
          uvT(loc,loc)+=delta;
          uvT/=opip;
          G=G*(Dmatrix::Identity(N2,N2)-uvT);
        }
      }
    }


    F=map(B+l*N4,N2,N2);
    G=F.colPivHouseholderQr().solve(G*F);


    //if lcount%some value==0 recalcG
    if(lcount%lcountmax==0){
      errG = buildG(G,B,L,N,m,(l+1)%L,lcount,avgerr,avgerrcounter);
      if(errG>10e-12&&lcountmax>=2) lcountmax--;
      if(errG<10e-12) lcountb12++;
      if(lcountb12==20){
        lcountb12=0;
        lcountmax++;
      }
    }
  }
}


int main(int argc, char *argv[]){

	//set up RNG
	mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);

	// read in parameters from command line
	double g,beta,w;
	int N,L,m=5,warmup,count,avgsteps,lagmax;
	int success=getparams(argc,argv,N,L,warmup,count,avgsteps,g,w,beta,lagmax);
	if(!success) return 1;
	double dT=beta/L;
	double lambda=g*g/8./w/w;
	double mu=-g*g/w/w;
	int N2=N*N;
	int N4=N*N*N*N;

	// read in the Local Neural Network Parameters
	Eigen::MatrixXd WL1(28,15);
	Eigen::MatrixXd WL2(15,5);
	Eigen::VectorXd WL3(5);
	double WLz;
	Eigen::VectorXd bL1(15);
	Eigen::VectorXd bL2(5);
	double bL3;
	Eigen::VectorXd xmean(28);
	Eigen::VectorXd xstd(28);
	double zmean;
	double zstd;
	
	
	std::ifstream in;
	in.open("./../L41local_hyper_parameter.dat");
	std::string word;
	cout<<in.is_open()<<endl;
	int i,j;
	for(i=0;i<2;i++) in>>word;
	for(i=0;i<28;i++){
		for(j=0;j<15;j++){
			in>>word;
			cout<<word;
			WL1(i,j)=atof(word.c_str());
		}
	}
	cout<<WL1<<endl;
	for(i=0;i<2;i++) in>>word;
	for(i=0;i<15;i++){
		for(j=0;j<5;j++){
			in>>word;
			WL2(i,j)=atof(word.c_str());
		}
	}
	for(i=0;i<2;i++) in>>word;
	for(j=0;j<5;j++){
		in>>word;
		WL3(j)=atof(word.c_str());
	}
	in>>word;
	in>>word;
	in>>word;
	WLz=atof(word.c_str());
	in>>word;
	for(j=0;j<15;j++){
		in>>word;
		bL1(j)=atof(word.c_str());
	}
	in>>word;
	for(j=0;j<5;j++){
		in>>word;
		bL2(j)=atof(word.c_str());
	}
	in>>word;
	in>>word;
	bL3=atof(word.c_str());
	for(j=0;j<28;j++){
		in>>word;
		xmean(j)=atof(word.c_str());
	}
	for(j=0;j<28;j++){
		in>>word;
		xstd(j)=atof(word.c_str());
	}
	in>>word;
	zmean=atof(word.c_str());
	in>>word;
	zstd=atof(word.c_str());
	in.close();

	//read in the global NN parameters
	double *Wconv=new double[3*8*10];
	Eigen::MatrixXd WG1(120,20);
	Eigen::VectorXd WG2(20);
	double WGg;
	double WGz;
	Eigen::VectorXd theta(10);
	Eigen::VectorXd bconv(10);
	Eigen::VectorXd bG1(20);
	double bG2;
	Eigen::MatrixXd xGmean(41,3);
	Eigen::MatrixXd xGstd(41,3);
	double dxGmean, dxGstd, zGmean, zGstd, gGmean, gGstd;

	in.open("L41global_hyper_parameter.dat");
	for(i=0;i<10;i++){
		for(j=0;j<8;j++){
			for(int k=0;k<3;k++){
				in>>word;
				Wconv[i*24+j*3+k]=atof(word.c_str());
			}
		}
	}
	for(i=0;i<120;i++){
		for(j=0;j<20;j++){
			in>>word;
			WG1(i,j)=atof(word.c_str());
		}
	}
	for(j=0;j<20;j++){
		in>>word;
		WG2(j)=atof(word.c_str());
	}
	in>>word;
	WGg=atof(word.c_str());
	in>>word;
	WGz=atof(word.c_str());
	for(j=0;j<10;j++){
		in>>word;
		theta(j)=atof(word.c_str());
	}
	for(j=0;j<10;j++){
		in>>word;
		bconv(j)=atof(word.c_str());
	}
	for(j=0;j<20;j++){
		in>>word;
		bG1(j)=atof(word.c_str());
	}
	in>>word;
	bG2=atof(word.c_str());
	for(i=0;i<41;i++){
		for(j=0;j<3;j++){
			in>>word;
			xGmean(i,j)=atof(word.c_str());
		}
	}
	for(i=0;i<41;i++){
		for(j=0;j<3;j++){
			in>>word;
			xGstd(i,j)=atof(word.c_str());
		}
	}
	in>>word;
	dxGmean=atof(word.c_str());
	in>>word;
	dxGstd=atof(word.c_str());
	in>>word;
	zGmean=atof(word.c_str());
	in>>word;
	zGstd=atof(word.c_str());
	in>>word;
	gGmean=atof(word.c_str());
	in>>word;
	gGstd=atof(word.c_str());
	in.close();


	//counters
	int lcount=0, lcountmax=5, lcountb12=0;
	int accept=0, NNtot=0, NNaccept=0;
	int total=0;
	int gtotal=0, gaccept=0;
	unsigned int avgerrcounter=0;

	//measurements
	double SCDW=0;
	double SCDW2=0;
	double SO=0;
	double avgerr=0;
	double *SCDWsamps = new double[lagmax];
	double *SCDWac = new double[lagmax];

	//make the e^K matrix
	//X is array of phonon positions first coordinate is time layer, second is y, third(continuous) is x
	//B are the B_l matricies there are L of them and each are N*NxN*N each stored rowmajor. within each row cont is x, discont is Y
	//G is the green's function it is an N*NxN*N matrix sorted in row-major manor.  within each row, the continuous index is X, discontinuous index is Y
	//buildG takes in G, all the B matricies, the number of B matricies, the number of sites, 
	//the number of products per UDR decomp, the l value that G is currently at and sets lcount to zero
	double *X = new double[L*N2];
	double *Xhold = new double[L*N2];
	double *B = new double[L*N4];
	double *Bhold = new double[L*N4];
	Eigen::MatrixXd G(N2,N2);


	{//fill the X and B and G matricies
	Dmatrix eK(N2,N2);
	Dmatrix K=Dmatrix::Zero(N2,N2);
	K.diagonal().array()-=mu;
	for(int y=0;y<N;y++){
		for(int x=0;x<N;x++){
			K(y*N+x,(y+1)%N*N+x)=-1;
			K(y*N+x,(y-1+N)%N*N+x)=-1;
			K(y*N+x,y*N+(x+1)%N)=-1;
			K(y*N+x,y*N+(x-1+N)%N)=-1;
		}
	}
	eK=(-dT*K).exp();
	for(int l=0;l<L*N2;l++){
		X[l]=unif(rng)*-2*g/w/w;
//		X[l]=-1;
	}
	Dmatrix Xmat=Dmatrix::Zero(N2,N2);
	int l,n;
	for(l=0;l<L;l++){
		for(n=0;n<N2;n++){
			Xmat(n,n)=-dT*g*X[l*N2+n];
		}
		map(B+(l*N4),N2,N2)=Xmat.exp()*eK;
	}
	}
	avgerr=0;
	avgerrcounter=0;

	double errG = buildG(G,B,L,N,m,0,lcount, avgerr, avgerrcounter);
	auto t1 = std::chrono::high_resolution_clock::now();

	cout<<"determinant"<<G.determinant()<<endl;

	//warmup updates
  for(int iter=0;iter<10000;iter++){
    //Local updates
    LocalUpdates(total, accept, errG, lcount, m, L, N, g, dT, w, unif, rng, X, G, B, lcountmax, lcountb12, avgerr, avgerrcounter);
    cout<<iter<<endl;
  }



	cout<<"determinant"<<G.determinant()<<endl;






	int largesugg=0;
	int largeacc=0;




	//warmup updates
	for(int iter=0;iter<warmup;iter++){
		double deaccum=0;
		memcpy(Xhold,X,sizeof(double)*N2*L);
		memcpy(Bhold,B,sizeof(double)*N4*L);
		//Local updates
		LocalUpdates(total, accept, errG, lcount, m, L, N, g, dT, w, unif, rng, X, G, B, lcountmax, lcountb12, avgerr, avgerrcounter, WL1, WL2, WL3, WLz, bL1, bL2, bL3, xmean, xstd, zmean, zstd,deaccum);

		//Global updates
	//	GlobalUpdates(gtotal, gaccept, count, lcount, m, L, N, g, dT, w, unif, rng, X, G, B, Wconv,	WG1, WG2, WGg, WGz, theta, bconv, bG1, bG2, xGmean, xGstd, dxGmean, dxGstd, zGmean, zGstd, gGmean, gGstd, deaccum);

		cout<<iter<<endl;
		double dS=calcS(X,w,dT,N2,L)-calcS(Xhold,w,dT,N2,L);
		errG = buildG(G,Bhold,L,N,m,0,lcount,avgerr,avgerrcounter);
		double numdet = G.determinant();
		errG = buildG(G,B,L,N,m,0,lcount, avgerr, avgerrcounter);
		double dendet = G.determinant();
		cout<<"deacc "<<deaccum<<endl;
		cout<<"dS"<<dS<<endl;
		cout<<"numdet"<<numdet<<endl;
		cout<<"dendet"<<dendet<<endl;
		cout<<"ratio"<<numdet*numdet/dendet/dendet<<endl;
		cout<<"The actual change in energy"<<exp(-dT*dS)*numdet*numdet/dendet/dendet<<endl;
		cout<<"The neural network's dE"<<exp(-beta*deaccum);
		double prob = exp(-dT*dS+beta*deaccum)*numdet*numdet/dendet/dendet;
		cout<<"prob"<<prob<<endl;
		if(unif(rng)>min(1.,prob)){
			memcpy(X,Xhold,sizeof(double)*N2*L);
			memcpy(B,Bhold,sizeof(double)*N4*L);
			largesugg++;
		}
		else{
			largesugg++;
			largeacc++;
		}
		cout<<(double)largeacc/largesugg;
	}
	total=accept=gtotal=gaccept=0;
	












	
	int x1,y1,x2,y2;
	double SCDWc, SOc;
	//averaging updates
	for(int iter=0;iter<avgsteps;iter++){

		double deaccum=0;
		memcpy(Xhold,X,sizeof(double)*N2*L);
		memcpy(Bhold,B,sizeof(double)*N4*L);

		//Local updates
		LocalUpdates(total, accept, errG, lcount, m, L, N, g, dT, w, unif, rng, X, G, B, lcountmax, lcountb12, avgerr, avgerrcounter, WL1, WL2, WL3, WLz, bL1, bL2, bL3, xmean, xstd, zmean, zstd, deaccum);
		cout<<iter<<endl;

		//Global updates
		GlobalUpdates(gtotal, gaccept, count, lcount, m, L, N, g, dT, w, unif, rng, X, G, B, Wconv,	WG1, WG2, WGg, WGz, theta, bconv, bG1, bG2, xGmean, xGstd, dxGmean, dxGstd, zGmean, zGstd, gGmean, gGstd, deaccum);

		
		cout<<iter<<endl;
		double dS=calcS(X,w,dT,N2,L)-calcS(Xhold,w,dT,N2,L);
		errG = buildG(G,Bhold,L,N,m,0,lcount,avgerr,avgerrcounter);
		double numdet = G.determinant();
		errG = buildG(G,B,L,N,m,0,lcount, avgerr, avgerrcounter);
		double dendet = G.determinant();
		double prob = exp(-dT*dS+beta*deaccum)*numdet*numdet/dendet/dendet;
		cout<<"prob"<<prob<<endl;
		if(unif(rng)>min(1.,prob)){
			memcpy(X,Xhold,sizeof(double)*N2*L);
			memcpy(B,Bhold,sizeof(double)*N4*L);
			largesugg++;
			errG=buildG(G,B,L,N,m,0,lcount,avgerr,avgerrcounter);
		}
		else{
			largesugg++;
			largeacc++;
		}
		cout<<(double)largeacc/largesugg;


		//Make Measurements
		SCDWc=0;
		SOc=0;
		for(y1=0;y1<N;y1++){
			for(x1=0;x1<N;x1++){
				SOc+=2*(1-G(y1*N+x1,y1*N+x1));
				SCDWc+=2.*G(y1*N+x1,y1*N+x1);
				for(y2=0;y2<N;y2++){
					for(x2=0;x2<N;x2++){
						SCDWc+=(-2*(abs(x1+y1-x2-y2)%2)+1)*(4*G(y1*N+x1,y1*N+x1)*G(y2*N+x2,y2*N+x2)-2*G(y1*N+x1,y2*N+x2)*G(y2*N+x2,y1*N+x1));
					}
				}
			}
		}
		SO=(iter*SO+SOc/N2)/(iter+1.);
		SCDW=(iter*SCDW+SCDWc/N2)/(iter+1.);
		SCDW2=(iter*SCDW2+SCDWc*SCDWc/N4)/(iter+1.);
		SCDWsamps[iter%lagmax]=SCDWc/N2;
		if(iter>=lagmax-1){
			for(int j=0;j<lagmax;j++){
				SCDWac[j]=((iter-lagmax+1.)*SCDWac[j]+SCDWsamps[(iter-(lagmax-1))%lagmax]*SCDWsamps[(iter-(lagmax-1)+j)%lagmax])/(iter-lagmax+2.);
			}
		}
	}

	auto t2 = std::chrono::high_resolution_clock::now();











	std::ofstream file;
	file.precision(15);
	file.open(to_string(L)+".dat");
	file<<"n: "<<N<<endl;
	file<<"g: "<<g<<endl;
	file<<"w: "<<w<<endl;
	file<<"beta: "<<beta<<endl;
	file<<"total time: "<<std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count()<<endl;
	file<<"warmup sweeps: "<<warmup<<endl;
	file<<"avg sweeps: "<<avgsteps<<endl;
	file<<"lcountmax: "<<lcountmax<<endl;
	file<<"averaged error: "<<avgerr<<endl;
	file<<"global acceptance prob "<<(double)accept/total<<endl;
	file<<"local acceptance prob "<<(double)gaccept/gtotal<<endl;
	file<<"large acceptance prob "<<(double)largeacc/largesugg<<endl;
	file<<"SO= "<<SO<<endl;
	file<<"SCDW= "<<SCDW<<endl;
	file<<"SCDW2= "<<SCDW2<<endl;
	file<<"LSCDW0= "<<SCDWac[0]<<endl;
	file<<"maxlag: "<<lagmax<<endl;
	for(int i=0;i<lagmax;i++){
		file<<(SCDWac[i]-SCDW*SCDW)/(SCDW2-SCDW*SCDW)<<endl;
	}
	
	delete[] SCDWac;
	delete[] X;
	delete[] SCDWsamps;
	delete[] B;
}
