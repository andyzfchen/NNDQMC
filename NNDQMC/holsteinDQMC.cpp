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
	//																QR decompose A
	//																set A to the orthogonal part, D to the diagonal of R, and R to rescaled R*old R
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


//This function performs the global update sweeps
void GlobalUpdates(int &gtotal, int &gaccept, int &count, int lcount, int &m, int &L, int &N, double &g, double &dt, double &w, std::uniform_real_distribution<double> &unif, mt19937_64 &rng, double *X, Dmatrix &G, double *B){

	int N2=N*N;
	int N4=N*N*N*N;

	//counters
	int i,l,k;
	double avge=0;
	unsigned int avgec=0;
	//values
	int sitex,sitey;
	Dmatrix deltaX(L,1);
	double dS,numdet,delta, *Bp, err, dendet;
	Dmatrix Gp(N2,N2);
	err=buildG(G,B,L,N,m,0,lcount,avge,avgec);

	//perform updates at count numer of sites
	for(i=0;i<count;i++){
		gtotal++;

		//the site we are doing the change on
		sitex=(int)floor(N*unif(rng));
		sitey=(int)floor(N*unif(rng));

		//suggest change and then calculate dS, detG
		for(l=0;l<L;l++){
			deltaX(l,0)=-2.*(X[l*N2+sitey*N+sitex])-2.*g/w/w;
			//dx=(unif(rng)>0.5)?-2*g/w/w:2*g/w/w;
		}
		dS=0;
		for(l=0;l<L;l++){
			dS+=2.*X[l*N2+sitey*N+sitex]+2.*g/w/w;
		}
		dS*=g;
		numdet=G.determinant();

		//change the B matricies
		for(l=0;l<L;l++){
			Bp=B+(l*N4+(sitey*N+sitex)*N2);
			delta=exp(-dt*g*deltaX(l,0))-1;
			for(k=0;k<N2;k++){
				Bp[k]*=(1+delta);
			}
		}

		//calculate new G from new Bs
		err = buildG(Gp,B,L,N,m,0,lcount,avge,avgec);
		dendet=Gp.determinant();

		//either accept or reject
		if(unif(rng)<exp(-dS*dt)*numdet*numdet/dendet/dendet){
			gaccept++;
			G=Gp;	
			for(l=0;l<L;l++){
				X[l*N2+sitey*N+sitex]+=deltaX(l,0);
			}
		}
		else{
			for(l=0;l<L;l++){
      				Bp=B+(l*N4+(sitey*N+sitex)*N2);
				delta=exp(-dt*g*deltaX(l,0))-1;
      				for(k=0;k<N2;k++){
        				Bp[k]/=(1+delta);
      				}
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


	//counters
	int lcount=0, lcountmax=5, lcountb12=0;
	int accept=0;
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
	double *B = new double[L*N4];
	Dmatrix G(N2,N2);

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
	double errG = buildG(G,B,L,N,m,0,lcount, avgerr, avgerrcounter);
	avgerr=0;
	avgerrcounter=0;

	auto t1 = std::chrono::high_resolution_clock::now();


	//warmup updates
	for(int iter=0;iter<warmup;iter++){
		//Local updates
		LocalUpdates(total, accept, errG, lcount, m, L, N, g, dT, w, unif, rng, X, G, B, lcountmax, lcountb12, avgerr, avgerrcounter);
		//Global updates
		GlobalUpdates(gtotal, gaccept, count, lcount, m, L, N, g, dT, w, unif, rng, X, G, B);
		cout<<iter<<endl;
	}
	total=accept=gtotal=gaccept=0;
	



	
	int x1,y1,x2,y2;
	double SCDWc, SOc;
	//averaging updates
	for(int iter=0;iter<avgsteps;iter++){
	//Local updates
		LocalUpdates(total, accept, errG, lcount, m, L, N, g, dT, w, unif, rng, X, G, B, lcountmax, lcountb12, avgerr, avgerrcounter);
		cout<<iter<<endl;
		//Global updates
		GlobalUpdates(gtotal, gaccept, count, lcount, m, L, N, g, dT, w, unif, rng, X, G, B);
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
	file<<"global acceptance prob"<<(double)accept/total<<endl;
	file<<"local acceptance prob"<<(double)gaccept/gtotal<<endl;
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
