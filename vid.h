#include <mpi.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

#define MAX_PROC 16
#define TAG 0
#define BUFSIZE 5

vector<int> process_heights(char **argv, int arg_index);
void debug_print_vector(vector<int> vect);
void debug_print_double_vector(vector<double> vect);
void debug_print_2D_vector(vector<vector<int>> vect);
bool cmdOptionExists(char** begin, char** end, const std::string& option);
vector<vector<int>> get_heights_per_proc(char **argv, int proc_cnt);
void send_heights(char **argv, vector<vector<int>> heights_per_proc);
int get_proc_cnt(char **argv);
double compute_vert_angle(int target, int base, int index);
int get_ht_id(int my_id, vector<vector<int>> heights_per_proc);
