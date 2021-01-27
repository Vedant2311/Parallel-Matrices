/*  The compile and run commands are:

	Compile:- $ g++-7 mr-pr-cpp.cpp /usr/lib/x86_64-linux-gnu/libboost_system.a /usr/lib/x86_64-linux-gnu/libboost_iostreams.a /usr/lib/x86_64-linux-gnu/libboost_filesystem.a -pthread -o mr-pr-cpp.o
	Run:- $ ./mr-pr-cpp.o <input_file>.txt -o <output_file>.txt
		
	* Compile command assumes the presence of the Mapreduce include files in the cwd
	* Run command assumes the existence of a file <input_file>.txt in a directory named "test" in the cwd of the code
	* Run command also assumes the existence of a directory named "output" where it writes the outputs of the pagerank algorithm execution.  

*/

// Including the Boost configurations
#include <boost/config.hpp>
#if defined(BOOST_MSVC)
#   pragma warning(disable: 4127)

// turn off checked iterators to avoid performance hit
#   if !defined(__SGI_STL_PORT)  &&  !defined(_DEBUG)
#       define _SECURE_SCL 0
#       define _HAS_ITERATOR_DEBUGGING 0
#   endif
#endif

#include <bits/stdc++.h>
#include "mapreduce.hpp"
using namespace std;

// Functions for the proper usage of the code
const char *OUTPUT_ARG = "-o";

void usage() {
    cerr << "Use the below format " << endl
    	 << "pagerank <graph_file> -o <output_file>" << endl
         << " -o enable output " << endl;
}

// Global variables for the pagerank algorithm
vector<vector<double>> intermediate_page(100000);
vector<vector<int>> outgoing_links(100000);
vector<double> pageranks(100000);
vector<double> pageranks_calc(100000);

// Parameters of alpha and convergence rate
double alpha = 0.85;
double conv = 0.00001;
double dangling_pointer =0.0f;
int total_webpages=0;

namespace page_rank {

template<typename MapTask>
class number_source : mapreduce::detail::noncopyable
{
  public:
    number_source()
      : sequence_(0)
    {
    }

    bool const setup_key(typename MapTask::key_type &key)
    {
        key = sequence_++;
        return (key <= total_webpages);
    }

    bool const get_data(typename MapTask::key_type const &key, typename MapTask::value_type &value)
    {
        value = pageranks[key];
        return true;
    }

  private:
    unsigned  sequence_;
};

struct map_task : public mapreduce::map_task<unsigned, double>
{
    template<typename Runtime>
    void operator()(Runtime &runtime, key_type const &key, value_type const &value) const
    {
		int chilren = outgoing_links[key].size();
		
		if(chilren!=0){
			for(int i =0; i<chilren; i++){
		        typename Runtime::reduce_task_type::key_type const emit_temp = outgoing_links[key][i];
		        double temp = double(value/chilren);
		        runtime.emit_intermediate(emit_temp, temp);			
			}
		}
		
		// Debugging
		runtime.emit_intermediate(key,0);
		
    }
};

struct reduce_task : public mapreduce::reduce_task<unsigned, double>
{
    template<typename Runtime, typename It>
    void operator()(Runtime &runtime, key_type const &key, It it, It ite) const
    {
        value_type results(*it);
		double pagerank_temp = 0.0f;
        for (It it1=it; it1!=ite; it1++)
        {
			pagerank_temp += (double)(*it1);
        }

		// Applying the pagerank update rule
		double pageranks_calc_temp = pagerank_temp*alpha + (1-alpha)/total_webpages + alpha*dangling_pointer/total_webpages;   
		
		// Emitting the (Key,Value) pair for the pageranks_calc
		results = pageranks_calc_temp;
		runtime.emit(key,results);
    }
};

typedef
mapreduce::job<page_rank::map_task,
               page_rank::reduce_task,
               mapreduce::null_combiner,
               page_rank::number_source<page_rank::map_task>
> job;

} // namespace page_rank


// Function to check the convergence
bool converging(vector<double> pageranks, vector<double> pageranks_calc, int total_webpages){
	bool isConverged= true;
	for (int i = 0; i < total_webpages; i++){
		if (abs(pageranks[i] - pageranks_calc[i])>conv){
			isConverged = false;
			break;
		}
	}
	return isConverged;
}

int main(int argc, char const *argv[])
{
	
	// File names
	string input_file;
	string output_file;
		
	// Reading the command line inputs		
	if (argc == 4){
		if (!strcmp(argv[2],OUTPUT_ARG)){
			input_file=argv[1];
			output_file=argv[3];
		}
		
		else{
			usage();
			exit(1);
		}
	}
	
	else{
		usage();
		exit(1);
	}
	
	int max_page_id = 0;
	int min_page_id = 0;	
	
	// Reading inputs from the specified file in the test folder of the standard pagerank repo
	ifstream fopen;
	fopen.open("test/" + input_file);
	while(!fopen.eof()){
		int a,b;
		fopen>>a>>b;
		outgoing_links[a].push_back(b);
		
		int max_temp = std::max(a,b);
		if (max_temp > max_page_id)
			max_page_id = max_temp;
		
		int min_temp = std::min(a,b);
		if (min_temp < min_page_id)
			min_page_id = min_temp;
	}
	
	// Total number of web-pages as: max_page_id - min_page_id + 1
	fopen.close();
	total_webpages = max_page_id - min_page_id + 1;

	// Initializing the pageranks;
	double pgr = double(1.0f/total_webpages);
	for(int i=0; i<total_webpages; i++){
		pageranks[i] = pgr;
	}

	// Initializing the Mapreduce specifications
    mapreduce::specification spec;
    spec.reduce_tasks = std::max(1U, std::thread::hardware_concurrency());

	while(true){
		
		cout << pageranks[0] << endl;
		
		// Handling dangling pointers
		dangling_pointer=0.0f;
		for(int i=0; i<total_webpages; i++){			
			if (outgoing_links[i].size()==0){
				dangling_pointer+=pageranks[i];
			}
		}		
		
		// Scheduling the map-reduce for the pagerank calculation
        page_rank::job::datasource_type number_source;
        page_rank::job job(number_source, spec);
        mapreduce::results result; 
        
		#ifdef _DEBUG
			job.run<mapreduce::schedule_policy::sequential<page_rank::job> >(result);
		#else
			job.run<mapreduce::schedule_policy::cpu_parallel<page_rank::job> >(result);
		#endif

		// Storing the reults of the calculated pagerank
        for (auto it=job.begin_results(); it!=job.end_results(); it++){
            pageranks_calc[it->first] = it->second;
        }    
		    
		// Checking for convergence
		if(converging(pageranks, pageranks_calc, total_webpages)){
			break;
		}	
		else{
			for(int i=0; i<total_webpages; i++){
				pageranks[i] = pageranks_calc[i];
			}		
		}
	}

	// Writing to the file
	ofstream fout;
	fout.open("output/" + output_file);

	string str;
	double sumVal = 0.0;
	for(int i=0; i<total_webpages; i++){
		fout<<i<<" = "<<pageranks[i]<<endl;
		sumVal += pageranks[i];
	}	
	
	fout << "sum " << sumVal << endl;
	fout.close();

	return 0;
}
