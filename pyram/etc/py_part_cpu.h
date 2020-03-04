extern void c_count_part(int* ndm_actual, int* nstar_actual,
                         int* nsink_actual, char* repository,
                         double* xmin, double* xmax,
                         double* ymin, double* ymax,
                         double* zmin, double* zmax,
                         int* cpu_list);



extern void c_load_part(double* star_float, int* star_int,
                        double*  dm_float, int* dm_int,
                        int* nstar_actual, int* ndm_actual,
                        int* nsink_actual, char* repository, 
                        double* xmin, double* xmax,
                        double* ymin, double* ymax,
                        double* zmin, double* zmax,
                        int* read_metal, int* cpu_list);

