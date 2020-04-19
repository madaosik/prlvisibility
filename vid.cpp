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
    MPI_Init(&argc,&argv);   
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);       
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);           

    //double start = MPI_Wtime();
    
    int total_in;
    vector<vector<int>> heights_per_proc = get_heights_per_proc(argv, numprocs, &total_in);

    int base = heights_per_proc[0][0];

    if (myid == 0) {
        //debug_print_2D_vector(heights_per_proc);
        send_heights(argv, heights_per_proc);
    }

    int my_heights_cnt = heights_per_proc[myid].size();
    MPI_Recv(mynums, my_heights_cnt, MPI_INT, 0, TAG, MPI_COMM_WORLD, &stat); 

    vector<double> angles(my_heights_cnt);
    
    if (myid == 0) {
        angles[0] = -2;
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
    double max_angle = find_max_in_array(angles, my_heights_cnt);

    // UPSWEEP faze
    int rounds= log2(numprocs);
    vector<double> values(rounds);
    int count = 0;
    int from, to, step, shift;
    double right;
    for (int round = 1; round <= rounds; round++)
    {
        values[round - 1] = max_angle;
        count++;
        step = pow(2, round);
        shift = pow(2, round - 1);
        if (myid % step == 0) {
            from = myid + shift;
            MPI_Recv(&right, 1, MPI_DOUBLE, from, TAG, MPI_COMM_WORLD, &stat);
            max_angle = max(max_angle, right);
            //cout << max_angle << endl;
        } else {
            to = myid - shift;
            MPI_Send(&max_angle, 1, MPI_DOUBLE, to, TAG, MPI_COMM_WORLD);
            break;
        }
    }

    // DOWNSWEEP faze
    double left;
    double new_value;
    int temp;

    max_angle = -2; // neutralni prvek
    
    for (int round = rounds; round > 0; round--)
    {
        shift = pow(2, round -1);
        step = pow(2, round);
        if (myid % step == 0) {
            to = myid + shift;
            if (myid == 0) {
                temp = 0;
            }
            else {
                temp = 1; 
            }
            left = values[--count - temp];
            new_value = max(left, max_angle);
            MPI_Send(&new_value, 1, MPI_DOUBLE, to, TAG, MPI_COMM_WORLD);
        }
        else if (round <= count) {
            from = myid - shift;
            MPI_Recv(&max_angle, 1, MPI_DOUBLE, from, TAG, MPI_COMM_WORLD, &stat);
        }
    }

    int visibility[my_heights_cnt];
    double new_max = max_angle;
    for (int i = 0; i < my_heights_cnt; i++)
    {   
        if (angles[i] > new_max) {
            new_max = angles[i];
            visibility[i] = true;
        }
        else {
            visibility[i] = false;
        }
    }

    // Finalni distribuce viditelnosti
    vector<int> result(total_in);
    int buf[BUFSIZE], sender_ht_cnt;
    for(int i = 1; i < numprocs; i++)
    {
        if (myid == i)
            MPI_Send(&visibility, my_heights_cnt, MPI_INT, 0, TAG, MPI_COMM_WORLD);
        if (myid == 0) {
            sender_ht_cnt = heights_per_proc[i].size();
            MPI_Recv(&buf, sender_ht_cnt, MPI_INT, i, TAG, MPI_COMM_WORLD, &stat);
            
            int offset = get_ht_id(i, heights_per_proc);
            for (int j = 0; j < sender_ht_cnt; j++) {
                result[j+offset] = buf[j];
            }
        }
    }

    if (myid == 0) {
        for (int i = 0; i < my_heights_cnt; i++)
            result[i] = visibility[i];
        //debug_print_vector(result);
        print_visibility(result);
    }
    /*
    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();
    if (myid == 0) {
        cout << end-start << endl;
    }
    */
    MPI_Finalize();
    return 0;
}

double find_max_in_array(vector<double> angles, int size)
{
    double max = angles[0];
    if (size == 1) return max;
    for (int i = 1; i < size; i++) {
        if (angles[i] > max) {
            max = angles[i];
        }
    }
    return max;
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

vector<vector<int>> get_heights_per_proc(char **argv, int proc_cnt, int* total_in)
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
    *total_in = pushed_cnt;
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
        if (log2(pow(2,i)) >= (inp_cnt / log2(pow(2,i))))
            break;

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
    cout << endl; 
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

void print_visibility(vector<int> vis)
{
    cout << "_";
    for (int i=1; i < vis.size(); i++) {
        if (vis[i] == 1) {
            cout << ",v";
        }
        else if (vis[i] == 0) {
            cout << ",u";
        }
    }
    cout << endl;
}