#import <string.h>
using namespace std;

map<string ,int> get_core_snapshot(string coredir, int snapshot, string sim){
	string fn = "m000p-"+snapshot;
	map<string, int> data;
	
	// ??? Alternative to h5py, how to read in files
	map<string, int> coredata = 
	
	for(int i=0; i<sizeof(keys)/sizeof(*keys); i++){
		data.insert(keys[i], coredata[keys[i]]);
	}

	return data;
}

map<string, int> add_coremass_column(map<string,int> corecat, string sim){
	// ??? How to mask
	map<string,int> mask = 
	return corecat
}

int get_parent_boundary(string[] foftags, string[] loc, int ncores, string name){

}

map<string, int> clean(map<string,int> corecat, map<string,int> sorted_coretags){

}

map<string, int> add_snapshot_to_trees(map<string,int> coretrees, map<string,int> corecat, int current_snap, string sim, int snap_index=0,int print_int=100, map<string,int> coretags_fofm=None, map<string,int> argsorted_fofm=None, map<string,int> sorted_indices=None,int ncore_min=0, int ncore_max=None, bool vector=False){

}

void write_binary(string outfile, map<string,int> coretrees, map<string,int> cores_to_write, map<string,int> foftags_to_write, bool vector=True, int start=None, int end=None, vector<int> column_counts=None){
}

int main(int argc, char* argv[]){

	// ??? C++ equivalent to argv.get()
	int nfiles = 
	int Nfiles = 
	int print_int =
	int ncore_min = 
	int ncore_max = 
	string sim =
	string name =
	string coredir = "../CoreCatalogs_"+sim;
	string treedir = "../CoreTrees/fof_group_"+sim;

	cout << "Simulation = " + sim << endl;
	cout << "Outputs written to " + treedir << endl;
	cout << "Writing " + nfiles " files out of a total of " Nfiles << endl;
	cout << "Reading core catalog from " + coredir << endl;
	
	// ??? How to read in filenames and parse
	string corefiles = 
	//int[] snapshots = 
	int[] snapshots = {}

	for(int i=0; i<sizeof(snapshots)/sizeof(*snapshots); i++){
		map<string, int> corecat = get_core_snapshot(coredir, snapshots[i], sim);
		if(corecat){
			corecat = add_coremass_column(corecat, sim);
			if(i==0){
				// ??? How to sort map by value
				map<string, int> sorted_coretags = 
				map<string, int> indices_fofm = 
				map<string, int> coretags_fofm = 
				map<string, int> foftags_fofm = 
				int Ncores = 
				int argsorted_coretags = 
				ncore_max = get_parent_boundary(
			}
			else{
				corecat = clean(corecat, sorted_coretags)
			}
		}
		coretrees = add_snapshot_to_trees(coretrees, corecat, int(s), sim, snap_index=n,coretrees = add_snapshot_to_trees(coretrees, corecat, int(s), sim, snap_index=n,coretags_fofm=coretags_fofm,argsorted_fofm=argsorted_coretags_fofm,sorted_indices=indices_fofm,print_int=print_int, ncore_max=ncore_max,ncore_min=ncore_min, vector=vector)
	}
	del corecat;

	int numcores = //length of keys
	int numfiles = nfiles;
	int stride = //numcores/numfiles
	int start = ncore_min;
	mode = ".serial";

	for(int i=0; i<nfiles; i++){
		string fn = treedir+string(i)+mode;
		string fn_bin = treedir+string(i)+mode;
		int end = min(start+stride, numcores);
		end = get_parent_boundary(foftags_fofm, end, Ncores, name="end");

		map<string,int> cores_to_write = // ??? How to slice
		map<string,int> foftags_to_write = // ??? How to slice

		if(sizeof(cores_to_write)/sizeof(*cores_to_write)>0){
			vector<int> column_counts;
			for(int i=0; i<sizeof(cores_to_write)/sizeof(*cores_to_write); i++)
				column_counts.push_back(/*len(coretrees[i]['Descendent']*/);
			write_binary(fn_bin, coretrees, cores_to_write, foftags_to_write, vector=vector, start=start, end=end, column_counts=column_counts);
		}
		start=end;
	}
	return;
}

