#include"nlib.h"
#include"math_nishi.h"
#include<Eigen/Dense>  // Eigen/Core, Geometry, LU, Cholesky, SVD, QR and Eigenvalues

using namespace Eigen;

/* ********************************************
 *    PCA from a file
 * ********************************************/
int pcanishi(  Inp_nishi inp1 ){

/* (0) variable definition
 */
   //int buf_int = 0;
   //unsigned int buf_unsigned_int = 0;
   //float buf_float = 0;
   double buf_double = 0;
   //long double buf_long_double = 0;

   vector<double> vec; //1D vector, including all input structural info. of all structures

/* (1) read input file
 * 
 */
   cout<<endl<<"--- INPUT INFORMATION --- \n";
   string infile = inp1.read("IHFILE"); // coord.dat   #this file should include one dimensional data for PCA
   //# output
   string outeigen = inp1.read("OUTEIGEN"); // out_eigen.txt
   string outpcmarker = inp1.read("OUTPCMARKER"); // _coord   #PC1, PC2, PC3, and so on will be output as named pc1"OUTPC".dat pc2"OUTPC".dat ...
   //# infomation
   unsigned int numofstru = atoi(inp1.read("NUMOFSTRU").c_str());
   unsigned int numofcomp = atoi(inp1.read("NUMOFCOMP").c_str());
   cout<<endl<<"Reading section of input file end. \n";

/* (2) read the file that include some structural informations
 *     from infile
 */
   cout<<endl<<"--- READING THE STRUCTURAL DATA FROM THE FILE --- \n";
   ifstream ifs( infile.c_str() ); // ifstream cannot handle type string (filename); it must be char*
   if(ifs.fail()){//error handling
      cerr<<"ERROR: Cannot open file: "<<infile<<endl;
      cerr<<"program was ended in reading file section of pcanishi.cpp\n";
      exit(1);
   }
   while(!ifs.eof()){
      ifs >> buf_double;
      vec.push_back(buf_double);
   }
   cout<<"vec[vec.size()-1] = "<<vec[vec.size()-1]<<endl;
   cout<<"Total num. of structures of read data  = "<<vec.size()/numofcomp<<endl;
   cout<<"Total num. of structures of your input = "<<numofstru<<endl;
   if(vec.size()/numofcomp != numofstru){
      cout<<"ERROR: read data and your input parameter are different\n";
      cout<<"program ended"<<endl;
      return 2;
   }

/* ********\\
 * 3-1   create vector Q = (q1x,q1y,q1z,q2x,...,q20z) 
 *     
 * */
   cout<<"\n--- PCA CALCULATION --- \n";
   //cout<<" 3-1  create vector Q \n";
   unsigned int dim_Q = numofcomp;
   cout<<"dimensionality of Q (components per one structure) is "<<dim_Q<<endl;
   unsigned int frame = numofstru;
   cout<<"Thr number of structures is "<<frame<<endl;

/* *********
 * 3-2   create variance-covariance matrix
 *
 * */
   cout<<" 3-2  create variance-covariance matrix \n";
   double *average_Q; double *average_QQ;
   average_Q = new double [dim_Q];
   average_QQ = new double [dim_Q*dim_Q];
   // initialize 2d-array
   for(unsigned int n=0;n<frame;n++){
      for(unsigned int i=0;i<dim_Q;i++){
         average_Q[i] = 0;
	for(unsigned int j=0;j<dim_Q;j++){
 	   average_QQ[ j + i * dim_Q ] = 0;  // 
 	}
      }
   }
   for(unsigned int n=0;n<frame;n++){
      for(unsigned int i=0;i<dim_Q;i++){
         average_Q[i] += vec[i+n*dim_Q];
			//cerr<<average_Q[i]<<" ";
			//cout<<"average_Q["<<i<<"] = "<<average_Q[i]<<endl;
	 for(unsigned int j=0;j<dim_Q;j++){
	    average_QQ[ j + i * dim_Q ] += vec[j+n*dim_Q] * vec[i+n*dim_Q];  // 
         }
      }
   }
        for(unsigned int i=0;i<dim_Q;i++){
                average_Q[i] = average_Q[i] / frame;
		//cerr<<"average_Q["<<i<<"] = "<<average_Q[i]<<endl;
        }
        for(unsigned int i=0;i<dim_Q*dim_Q;i++){
                average_QQ[i] = average_QQ[i] / frame;
		//cout<<"average_QQ["<<i<<"] = "<<average_QQ[i]<<endl;
        }
	//delete[] Q;
	
        MatrixXd VCV(dim_Q,dim_Q);// Eigen/Dense
	VCV = MatrixXd::Zero(dim_Q,dim_Q); //initialize matrix VCV with 0
        //Matrix<double, dim_Q,dim_Q> VCV;//not work
        for(unsigned int n=0;n<dim_Q;n++){
                for(unsigned int i=0;i<dim_Q;i++){
                	VCV(i,n) = average_QQ[ i + n * dim_Q ] - average_Q[i] * average_Q[n];
                	//cout<<"\nVCV(i,n) = "<<average_QQ[ i + n * dim_Q ]<<" - "<<average_Q[i] * average_Q[n];
			//cout<<" "<<VCV(i,n);
			//cerr<<VCV(i,n)<<" ";
		}
        	//cout<<endl;
	}
	//delete average_Q;
	delete[] average_QQ;
	//cerr<<VCV<<endl;
        for(unsigned int n=0;n<dim_Q;n++){
                for(unsigned int i=0;i<dim_Q;i++){
			//cout<<"DEBUG: VCV("<<i<<","<<n<<") = "<<VCV(i,n)<<", VCV("<<n<<","<<i<<") = "<<VCV(n,i)<<endl;
			if( VCV(i,n) != VCV(n,i) )cout<<"ERROR: !!!!!!!!!!!! VCV is not symmetric !!!!!!!!!!!!!!\n";
                }
        }
/* **********
 * 3-3   calcuate eigen-value and eigen-vector
 *
 * */

/*   Matrix2f A; // Eigen2 test section
   A << 1,2,2,3;
   SelfAdjointEigenSolver<Matrix2f> es1( A );
   cout<<"Eigen value: \n"
       << es1.eigenvalues() << endl;
   cout<<"Eigen vector: \n"
       << es1.eigenvectors() << endl;
*/
   cout<<" 3-3  calculate eigen-value \n";
	//SelfAdjointEigenSolver< MatrixXd > es( A );
	SelfAdjointEigenSolver< MatrixXd > es( VCV );
        if (es.info() != Success) abort();
	cout<<"eigenvalues in ascending order\n"
	    <<es.eigenvalues()<<endl<<endl;
	//cout<<"!!!!\n"<<es.eigenvalues().rows()<<", dim_Q = "<<dim_Q<<endl;
	//cout<<"\neigenvectors\n"<<es.eigenvectors()<<endl<<endl;
	//cout<<"The maximum eigenvalue of VCV\n"<<es.eigenvalues()[dim_Q - 1]<<endl;
	//cout<<"The maximum eigenvalue of VCV\n"<<es.eigenvalues().row(dim_Q - 1)<<endl;
	//cout<<"The eigenvector of the maximum eigenvalue\n"<<es.eigenvectors().col(dim_Q-1).transpose()<<endl;
	//cout<<"The second eigenvalue of VCV\n"<<es.eigenvalues()[dim_Q-2]<<endl;
	//cout<<"The eigenvector of the second eigenvalue\n"<<es.eigenvectors().col(dim_Q-2).transpose()<<endl;
   //cout<<"DRBUG: mark"<<endl;
  
   ofstream ofs_eigen;
   ofs_eigen.open( outeigen.c_str() );
   ofs_eigen<<"> Eigen values\n"<<es.eigenvalues()<<endl;
 
   for(unsigned int n=0;n<dim_Q;n++){
      ofs_eigen<<"> Eigen vector no"<<n+1<<"\n"<<es.eigenvectors().col(n).transpose()<<endl;
   }
   
   ofs_eigen<<"> average Q"<<endl;
   for(unsigned int i=0;i<dim_Q;i++){
      ofs_eigen<<average_Q[i]<<endl;
   }

   ofs_eigen.close(); 
   cout<<"ouput "<<outeigen<<"\n";

/* **********
 *  3-5  calculate the coordinates of the ensemble along the PCA axes
 *
 * */
   cout<<" 3-5  calculate components along principle axis \n";

	//double c1[frame][dim_Q]; 

	vector<double> buf_vec;  //for mapping, average_Q -> bef_vec -> Vector Q2
        for(unsigned int i=0;i<dim_Q;i++){
	        buf_vec.push_back(average_Q[i]);
        }
	VectorXd Q2=Map<VectorXd>(&buf_vec[0],buf_vec.size());
	//delete[] average_Q;
	//cout<<"Vector Q2 was set"<<endl;

   double c1[frame]; 
   for(unsigned int j=0;j<dim_Q;j++){
        for(unsigned int n=0;n<frame;n++){  
		buf_vec.clear();
                for(unsigned int i=0;i<dim_Q;i++){
			buf_vec.push_back(vec[i+n*dim_Q]);
		}
		VectorXd Q1 = Map<VectorXd>(&buf_vec[0],buf_vec.size());  //buf_vec -> Q1, raw data
		//cout<<"Vector Q1 was set"<<endl;
	        //for(unsigned int i=0;i<dim_Q;i++){ // Q - <Q>
        	c1[n] = es.eigenvectors().col(dim_Q - j - 1).transpose()*( Q1 - Q2 );
        	//cout<<"component = "<<c1[n]<<" , "<<n+1<<endl;
	}	
	//cout<<"PC was set"<<endl;

/* (4) output OUT
*/
   //for(unsigned int i=0;i<dim_Q;i++){
      string pc_num;  char buf[256];
      sprintf(buf,"PC%d%s.dat",j+1,outpcmarker.c_str());  //itoa()
      pc_num = buf;
      cout<<"\noutput section of "<<pc_num<<endl;

	ofstream ofs;
		ofs.open( pc_num.c_str() );
	        for(unsigned int n=0;n<frame;n++){
        		ofs<<c1[n]<<endl;
		}
	ofs.close();
   	cout<<"output \n";
	ofs.close();
   }
/* (6) calculate contribution ratio of c1 and c2
 *     
 * */
   cout<<endl<<"REPORT> (6) calculate contribution ratio \n";
   double sum_lambda = 0;
   for(unsigned int i=0;i<dim_Q;i++){
      sum_lambda += es.eigenvalues()[i];
   }
   cout<<"Summation of eigen values = "<<sum_lambda<<endl<<endl;
   double count_cntr_rt = 0;
   for(unsigned int i=0;i<dim_Q;i++){ 
      double cntr_rt;  //CoNTRibution_RaTion
      cntr_rt = es.eigenvalues()[dim_Q - 1 -i] / sum_lambda; 
      count_cntr_rt += cntr_rt;
      cout<<"contribution ratio (PC"<<i+1<<") = "<<cntr_rt<<endl;
   }
   cout<<"\nSummation of cntribution ratio = "<<count_cntr_rt<<endl;

// END
        return 0;
}
