//#pragma omp parallel for // parallel execution by OpenMP
for (int mu = 3; mu < 5; mu++) { //mu_f_max = 20
    mu_f = (double) mu;
    for (int sigm = 3; sigm < 4 ; sigm++) { //sigm_f_max = 0.7
        sigm_f = (1.0 / (double) sigm) * 0.9;
        special_face_design = DCC_Processing(special_faces_sequence, State_sVector, P_type, mu_f, sigm_f, RW_series_vector);
    } //  for (int sigm = 1; sigm < 10 ; sigm++)
} // for (int mu = 2; mu < 3; mu++)

/// ===== Elapsing time Processing ================
//        if (ProcessingON(confpath, time_step_one)) {
unsigned int Processing_time = clock();
P_time = (double) Processing_time;
cout << "Processing time is equal to  " << P_time / pow(10.0, 6.0) << "  seconds" << endl; cout << "-------------------------------------------------------------------------" << endl;
Out_logfile_stream << "Processing time is equal to  " << P_time / pow(10.0, 6.0) << "  seconds" << endl; Out_logfile_stream << "-------------------------------------------------------------------------" << endl;