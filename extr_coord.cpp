#include"nlib.h"
#include"math_nishi.h"
#include<Eigen/Dense>  // Eigen/Core, Geometry, LU, Cholesky, SVD, QR and Eigenvalues
//#include"inpnishi.h"

using namespace Eigen;

// constant values
#define GAS_CONST 8.31451  // gas constant, R ( joule/mol*k )
//#define BOLTZMAN_CONST 1.380658e-23  // Boltzman constant, kb ( joule/k )
#define BOLTZMAN_CONST 8.3144621  // Boltzman constant, kb ( joule/k*mol )
// kb = R / Na, where Na is Avogadro constant (6.02214129e-23)
// printf("Boltzman constant kb = %e \n", kb);
#define JOULE_CALORIE 4.184 // joule-calorie conversion unit

#define DEBUG_PCA 1 // for debugging, instead of comment out
//      #if DEBUG_PCA == 1
//            ........................
//      #endif


vector<double> quaternion( vector<double> &vec_ref, vector<double> &vec_tar );  // give two vectors as pointer, so please notice the changes of these vectors in this function will remain out of it.
int transfer_quat( vector<double> &vec, vector<double> &transf );
int rotate_quat( vector<double> &vec, vector<double> &rot_mat );


int main(int argc, char *argv[]){
  cout<<"Version info. extr_coord v1.0.0 \n";
  cout<<"Description: extract x,y,z coordinates of selected residues from trajectry file\n";
// ##################### ARGUMENT HANDLING ##########################
// argv[1]: input parameter file
  if( argv[1]==NULL ){
    puts("No ARGUMEMTS");
    puts("USAGE: ./extr_coord.exe (argv[1]: input parameter file)" );
    return 1;
  }
  cout<<"Your input-parameter file: "<<argv[1]<<endl;

// INPUT_PARAMETERS
   Inp_nishi inp1( argv[1] );
   
// ############# trajectory  ########################################
/* (1) load trajectory by tra_nishi
 * 
 */
   //cout<<endl<<"REPORT> (1) load trajectory \n";
   cout<<endl<<"--- TRAJECTORY INFORMATION --- \n";
   string codname = inp1.read("COD1");
   string pdbname = inp1.read("REFPDBNAME");
   string pcaatom = inp1.read("PCAATOM") ;
   int stride = atoi(inp1.read("STRIDE").c_str());
   int startframe = atoi( inp1.read("STARTFRAME").c_str() ) - 1 ;
   string superpbase = inp1.read("SUPERPBASE");   //v.1.1.0

   if(stride <= 0){
      return -1;
   }
   cout<<"loading\n"<<codname<<" and "<<pdbname<<endl;
	tra_nishi* tra1;
	tra1 = new tra_nishi(codname.c_str(), pdbname.c_str(), stride, pcaatom);
	cout<<"TOTAL FRAME = "<<tra1->total_step<<endl;
	cout<<"TOTAL ATOM = "<<tra1->pdb1->total_atom<<endl;
	cout<<"TOTAL SELECTED ATOM = "<<tra1->total_sel<<endl;
	//unsigned int frame = tra1->total_step;
/* (2) setting of region
 *     from intra_start to intra_end (internal number of atom)
 */
   //cout<<endl<<"REPORT> (2) specify the region \n";
   cout<<endl<<"--- RESIDUE RANGE --- \n";
   char startchain = inp1.read("STARTCHAIN").c_str()[0];
   char endchain = inp1.read("ENDCHAIN").c_str()[0];
   int startres = atoi(inp1.read("STARTRES").c_str());
   int endres = atoi(inp1.read("ENDRES").c_str());
   int intra_start, intra_end;
	intra_start = tra1->pdb1->search_n( startchain , startres );
	intra_end = tra1->pdb1->search_n_end( endchain , endres );
	//cout<<"intra_start = "<<intra_start<<endl;
	//cout<<"intra_end = "<<intra_end<<endl;
	if( intra_start < 0 || intra_end < 0 ){
           cout<<"ERROR: selection is not correct\n";
	   return -1;
	}

	cout<<"unselected: ";tra1->pdb1->disp_line(intra_start-1);
	cout<<"  selected: ";tra1->pdb1->disp_line(intra_start);
	cout<<"  selected: ";cout<<"...\n";
        cout<<"  selected: ";tra1->pdb1->disp_line(intra_end);
        cout<<"unselected: ";tra1->pdb1->disp_line(intra_end+1);

        //tra1->pdb1->write_pdb("zzz.pdb"); cout<<"output pdb file (zzz.pdb)\n";

/* (3) ATOM SELECTION
*/
   cout<<"\n--- ATOM SELECTION --- \n";
   int rej_ca=0, rej_h=0, rej_wat=0, rej_cim=0, rej_cip=0, rej_mainchain=0;
   int rtrn_sel, count=0, start_sel = 0, end_sel = 0;
   vector<double> vec_ref;
      for( int i=0;(signed)i<(signed)tra1->pdb1->total_atom;i++){
      //for( int i=intra_start;i<=intra_end;i++){
         rtrn_sel = select_atom( *tra1->pdb1, vec_ref, pcaatom, i );
	 if( i>=intra_start && i<=intra_end ){
	    switch( rtrn_sel ){
	    case 0: 
	       if( start_sel == 0){
	          start_sel = count;
		  tra1->pdb1->disp_line(i);
   	       }
	       end_sel = count ;
	       tra1->pdb1->disp_line(i);
	       break;
	    case 1: rej_ca++; break;
	    case 2: rej_mainchain++; break;
	    case 3: rej_h++; break;
	    case 4: rej_wat++; break;
	    case 5: rej_cim++; break;
	    case 6: rej_cip++; break;
	    default: cout<<"Unknown value of rtrn_sel \n";
	    }
	 }
	 if(rtrn_sel==0)count++;
      }
      
      //int cycle_frame = end_sel - start_sel + 1;
      //cout<<"!!!! start_sel and end_sel = "<<start_sel<<" and "<<end_sel<<endl;
      cout<<"Atoms in region between STARTRES and ENDRES = "<<intra_end - intra_start + 1<<endl;
      cout<<"Num. of atoms selected for PCA calculation = "<<end_sel - start_sel +1<<" in reference"<<endl;
      cout<<"Rejected atoms are as follows (Selected: "<<pcaatom<<")\n";
      cout<<"!CA atoms = "<<rej_ca<<endl;
      cout<<"!CA & !N & !C & !O atoms = "<<rej_mainchain<<endl;
      cout<<"element H (Hydrogens) = "<<rej_h<<endl;
      cout<<"residue name WAT (Water molecules) = "<<rej_wat<<endl;
      cout<<"residue name CIM (minus ions) = "<<rej_cim<<endl;
      cout<<"residue name CIP (plus ions) = "<<rej_cip<<endl;
   
   //vec.clear();
   int cod_num_i=2;
   cout<<"\n--- LOADING STRUCTURE ENSEMBLES --- \n";
   cout<<"Reading COD1"<<endl;

   vector<double> vec, vec_buf;
   if( superpbase == "YES" ){
      for(int i=0;i<start_sel*3;i++){
        vec_buf.push_back(vec_ref[i]);
      }
      for(unsigned int i=end_sel*3+3;i<vec_ref.size();i++){
        vec_buf.push_back(vec_ref[i]);
      }
   }
   else{
      for(int i=start_sel*3;i<=end_sel*3+2;i++){
        vec_buf.push_back(vec_ref[i]);
      }
   }
   //cout<<"!!!!!!! vec_buf.size() = "<<vec_buf.size()<<endl;
   //cout<<"!!!!!!! vec_ref.size() = "<<vec_ref.size()<<endl;
   //cout<<"!!! vec_ref[0] = "<<vec_buf[0]<<", vec_ref[n] = "<<vec_buf[vec_buf.size() -1]<<endl;

   string flag_rotmat = "notyet"; vector<double> rot_mat_1st;  //v.1.1.2
flag100:
   cout<<"TOTAL FRAME = "<<tra1->total_step<<endl;
   cout<<"TOTAL ATOM = "<<tra1->pdb1->total_atom<<endl;
   cout<<"TOTAL SELECTED ATOM = "<<tra1->total_sel<<endl;
   vector<double> vec_tar, rot_mat, transf;
   //cout<<"!!!!! tra1->cordx.size() = "<<tra1->cordx.size()<<endl;
   for(unsigned int n=startframe;n<tra1->total_step;n++){
      vec_tar.clear();
      if( superpbase == "YES" ){
         for(int i=0;i<start_sel;i++){
            vec_tar.push_back(tra1->cordx[n*tra1->total_sel+i]);
            vec_tar.push_back(tra1->cordy[n*tra1->total_sel+i]);
            vec_tar.push_back(tra1->cordz[n*tra1->total_sel+i]);
         }
         for(unsigned int i=end_sel+1;i<tra1->total_sel;i++){
            vec_tar.push_back(tra1->cordx[n*tra1->total_sel+i]);
            vec_tar.push_back(tra1->cordy[n*tra1->total_sel+i]);
            vec_tar.push_back(tra1->cordz[n*tra1->total_sel+i]);
         }
         //cout<<"!!!!!!! vec_tar.size() = "<<vec_tar.size()<<endl;
         //cout<<"!!!!!!! total_sel*3 = "<<tra1->total_sel*3<<endl;
         rot_mat = quaternion( vec_buf, vec_tar );  // give two vectors as pointer, so please notice the changes of these vectors in this function will remain out of it.
	 if(flag_rotmat == "notyet"){
	    rot_mat_1st = rot_mat;
	    flag_rotmat = "alreadydone";
	 }
	 vec_tar.clear();
         for(int i=start_sel;i<=end_sel;i++){
            vec_tar.push_back(tra1->cordx[n*tra1->total_sel+i]);
            vec_tar.push_back(tra1->cordy[n*tra1->total_sel+i]);
            vec_tar.push_back(tra1->cordz[n*tra1->total_sel+i]);
         }
         transf.clear();  // transfer target
         transf.push_back( rot_mat[12] );
         transf.push_back( rot_mat[13] );
         transf.push_back( rot_mat[14] );
         transfer_quat( vec_tar, transf );

         rotate_quat( vec_tar, rot_mat );  // rotate target
       
      }
      else {
         for(int i=start_sel;i<=end_sel;i++){
            vec_tar.push_back(tra1->cordx[n*tra1->total_sel+i]);
            vec_tar.push_back(tra1->cordy[n*tra1->total_sel+i]);
            vec_tar.push_back(tra1->cordz[n*tra1->total_sel+i]);
         }
         ///cout<<n<<" rmsd before = "<<rmsd2(vec_tar,vec_buf)<<endl;
         rot_mat = quaternion( vec_buf, vec_tar );  // give two vectors as pointer, so please notice the changes of these vectors in this function will remain out of it.
	 if(flag_rotmat == "notyet"){
	    rot_mat_1st = rot_mat;
	    flag_rotmat = "alreadydone";
	 }
      }
      //cout<<"!!! vec_tar[0] = "<<vec_tar[0]<<", vec_tar[n] = "<<vec_tar[vec_tar.size() -1]<<endl;
      //cout<<n<<" rmsd after = "<<rmsd2(vec_tar,vec_buf)<<endl;
      //cout<<"!!! vec_ref[0] = "<<vec_buf[0]<<", vec_ref[n] = "<<vec_buf[vec_buf.size() -1]<<endl;
      //cout<<"!!! vec_tar[0] = "<<vec_tar[0]<<", vec_tar[n] = "<<vec_tar[vec_tar.size() -1]<<endl;
   
      //cout<<"!!!!! size of cycle_frame and vec_tar = "<<cycle_frame<<" and "<<vec_tar.size()<<endl;
      for(unsigned int i=0;i<vec_tar.size();i++){
         vec.push_back( vec_tar[i] );
      }
   }
   //cout<<"!!!!! size of vec_ref and vec_tar = "<<vec_buf.size()<<" and "<<vec_tar.size()<<endl;
   //tra1->write_step("zzz.pdb",tra1->total_step -1);

   string cod_num;  char buf[32];
   sprintf(buf,"COD%d",cod_num_i);  //itoa()
   cod_num = buf;
   cout<<"\nReading "<<cod_num<<endl;
   cod_num_i++;
   codname = inp1.read(cod_num);
   if( codname == "nothing" ){
      cout<<"Section of reading trajectories was ended normaly\n";
      delete tra1;
   }else{
      delete tra1;
      tra1 = new tra_nishi(codname.c_str(), pdbname.c_str(), stride, pcaatom);
      goto flag100;
   }
   //unsigned int dim_0 = end_sel - start_sel +1;

   cout<<"\n--- READING REFERENCE STRUCTURE --- \n";
/*   cout<<"rot_mat[9] = "<<rot_mat[9]<<endl;
   cout<<"rot_mat[10] = "<<rot_mat[10]<<endl;
   cout<<"rot_mat[11] = "<<rot_mat[11]<<endl;
   cout<<"rot_mat_1st[9] = "<<rot_mat[9]<<endl;
   cout<<"rot_mat_1st[10] = "<<rot_mat[10]<<endl;
   cout<<"rot_mat_1st[11] = "<<rot_mat[11]<<endl;*/
   transf.clear();  //no need to transfer referencer because do not calculate RMSD
         transf.push_back( rot_mat_1st[ 9] );
         transf.push_back( rot_mat_1st[10] );
         transf.push_back( rot_mat_1st[11] );
         transfer_quat( vec_ref, transf );
	 for(int i=start_sel*3;i<=end_sel*3+2;i++){
	    vec.push_back(vec_ref[i]);
	 }
	 //cout<<"!!! vec.size() = "<<vec.size()<<endl;
	 //cout<<"!!! vec_ref.size()/3 = "<<vec_ref.size()/3<<endl;
   cout<<"The final Data of the structure set is of the Reference. \n"<<endl;

   cout<<"\n--- OUTPUT SECTION --- \n";
   //cout<<endl<<"REPORT> (3) PCA calculation starts \n";
/* ********\\
 * 3-1   create vector Q = (q1x,q1y,q1z,q2x,...,q20z) 
 *     
 * */
   //cout<<" 3-1  create vector Q \n";
   unsigned int dim_Q = vec_tar.size();
   cout<<"dimensionality of Q (components per one structure) is "<<dim_Q<<endl;
   //cout<<"cycle_frame = "<<cycle_frame<<endl; //cycle_frame should be dim_Q/3
   unsigned int frame = vec.size() / dim_Q;
   cout<<"Thr number of structures is "<<frame<<endl;

/* *********
 * output section
 *
 * */
   
   string outfile = inp1.read("OUTFILE");
   ofstream ofs;
   ofs.open( outfile.c_str() );
   for(int i=0;i<frame;i++){
      for(int j=0;j<dim_Q;j++){
         ofs<<vec[i*dim_Q + j]<<"    ";
      }
      ofs<<endl;
   }
   ofs.close();

// END
        return 0;
}
