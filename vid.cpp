/*
* Viditelnost - PRL projekt 3
* Author:      Adam Lanicek
*/

#include "vid.h"

int main(int argc, char *argv[])
{
    int numprocs;               // pocet procesoru
    int myid;                  
    MPI_Status stat;            
    int mynums[BUFSIZE];

    /* Program spusten pro ucely zjisteni poctu potrebnych procesoru */
    if(cmdOptionExists(argv, argv+argc, "-gp"))
    {
        cout << get_proc_cnt(argv) << endl;
        return 0;
    }

    /* Inicializace MPI */
    MPI_Init(&argc,&argv);      // inicializace MPI 
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);       // ulozime si celkovy pocet rozbehnutych procesu
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);           // id procesu

    vector<vector<int>> heights_per_proc = get_heights_per_proc(argv, numprocs);
    int base = heights_per_proc[0][0];

    if (myid == 0) {
        //debug_print_2D_vector(heights_per_proc);
        send_heights(argv, heights_per_proc);
    }

    int my_heights_cnt = heights_per_proc[myid].size();
    MPI_Recv(mynums, my_heights_cnt, MPI_INT, 0, TAG, MPI_COMM_WORLD, &stat); //buffer,velikost,typ,rank odesilatele,tag, skupina, stat

    vector<double> angles(my_heights_cnt);
    
    if (myid == 0) {
        angles[0] = 0;
        for (int i = 1; i < my_heights_cnt; i++) 
        {
         angles[i] = compute_vert_angle(mynums[i], base, i);
        }
    } else {
        int base_iter = get_ht_id(myid, heights_per_proc);
        int ht_iter = base_iter;
        for (int i = 0; ht_iter < (base_iter + my_heights_cnt); i++, ht_iter++)
        {
         angles[i] = compute_vert_angle(mynums[i], base, ht_iter);
        }
    }
    /*
    if (myid == 0) {
        debug_print_double_vector(angles);
    }
    if (myid == 3) {
        debug_print_double_vector(angles);
    }
    */




    //cout << numprocs << endl << endl;
    //double start, end;
    //MPI_Barrier(MPI_COMM_WORLD);
    //start = MPI_Wtime();
    MPI_Finalize(); 
    return 0;
}

int get_ht_id(int my_id, vector<vector<int>> heights_per_proc)
{
    if (my_id == 0) return 1;
    int ht_index = 0;
    for (int i=0; i < my_id; i++)
    {
        ht_index += heights_per_proc[i].size();
    }
    return ht_index;
}

double compute_vert_angle(int target, int base, int index)
{
    return atan((target - base)/index);
}

vector<vector<int>> get_heights_per_proc(char **argv, int proc_cnt)
{
    vector<vector<int>> heights_proc(proc_cnt);
    vector<int> heights = process_heights(argv,1);
    int heights_cnt = heights.size();
    int heights_per_cpu = heights_cnt / proc_cnt;
    int rem_heights = heights_cnt % proc_cnt; // zbyva prerozdelit, max hodnota proc_cnt - 1

    int height_iter = 0;
    int pushed_cnt = 0;
    int iter_cnt;
    vector<int> temp_heights;
    for (int i = 0; i < proc_cnt; i++)
    {
        iter_cnt = 0;
        for (int j = 0; j < heights_per_cpu; j++)
        {   
           temp_heights.push_back(heights[pushed_cnt + j]);
           iter_cnt++;
        }
        if (i < rem_heights) {
            temp_heights.push_back(heights[pushed_cnt + heights_per_cpu]);
            iter_cnt++;
        }
        heights_proc[i] = temp_heights;
        temp_heights.clear();
        pushed_cnt += iter_cnt;
    }
    return heights_proc;
}

void send_heights(char **argv, vector<vector<int>> heights_per_proc)
{
    int proc_id = 0;
    int len = heights_per_proc.size();
    int int_buf[BUFSIZE];

    for (int i = 0; i < len; i++) {
        int inner_len = heights_per_proc[i].size();
        for (int j = 0; j < inner_len; j++)
        {
            int_buf[j] = heights_per_proc[i][j];
        }
        MPI_Send(int_buf, heights_per_proc[i].size(), MPI_INT, i, TAG, MPI_COMM_WORLD);
    }
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

int get_proc_cnt(char **argv)
{
    vector<int> heights = process_heights(argv,2);
    int inp_cnt = heights.size();
    int proc_cnt = 1;

    int i = 1;
    while (log2(proc_cnt) < (inp_cnt / proc_cnt))
    {   
        proc_cnt = pow(2, i);
        if (proc_cnt == MAX_PROC)
            break;
        i++;
    }
    return proc_cnt;
}

vector<int> process_heights(char **argv, int arg_index)
{
    string str = argv[arg_index];
    vector<int> heights;

    stringstream ss(str);

    for (int i; ss >> i;) {
        heights.push_back(i);    
        if (ss.peek() == ',')
            ss.ignore();
    }
    return heights;
}

void debug_print_vector(vector<int> vect)
{
    for (std::size_t i = 0; i < vect.size(); i++)
        cout << vect[i] << std::endl;
    cout << endl; 
}

void debug_print_double_vector(vector<double> vect)
{
    for (std::size_t i = 0; i < vect.size(); i++)
        cout << vect[i] << std::endl;
    //cout << endl; 
}

void debug_print_2D_vector(vector<vector<int>> vect)
{
    int firstD_size = vect.size();
    int secondD_size;

    for (int i = 0; i < firstD_size; ++i)
    {
        secondD_size = vect[i].size();
        for (int j = 0; j < secondD_size; ++j)
        {
            cout << vect[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;
}