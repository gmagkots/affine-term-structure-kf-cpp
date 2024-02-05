void vasicek::write_hdf5()
{
    hid_t   file_id, group_id, dataset_id;
    hid_t   dataspace_id, dataspace_id2, memspace_id;
    hsize_t dims[1], dims2[1];
    herr_t  status;
    stringstream int2str;
    string groupname;
    double *h5feeder;
    double scratch[1];

    dims2[0] = 1;

    // Create a new file using default properties
    file_id = H5Fcreate("bond_output.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Create the dataspace for the scalar T
    dataspace_id2 = H5Screate_simple(1, dims2, NULL);

    for (size_t bonds=0; bonds<get_nbond(); bonds++){

        /* Uncustommary: Create multiple dataspaces (each for a bond),
           because each bond contains different amount of quoted prices. */
        dims[0] = zero_rate[bonds].size();
        //dims[1] = 1; 
        dataspace_id = H5Screate_simple(1, dims, NULL);

        /* Allocate the native array which will take data
           from array2D's and feed it to the hdf5 routines. */
        h5feeder = new double [dims[0]];

        // Create the group corrsponding to a bond
        int2str << bonds; groupname = "Bond";
        groupname += int2str.str();
        group_id = H5Gcreate(file_id, groupname.c_str(), H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);

        // Create a dataset to include time to maturity T
        /*scratch[0] = bond_maturity_T[bonds];
        dataset_id = H5Dcreate(group_id, "time_to_maturity_T", H5T_NATIVE_DOUBLE,
                                 dataspace_id2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          scratch);
                          status = H5Dclose(dataset_id);*/

        // Write the data
        memcpy( h5feeder, &h5_time[bonds][0], sizeof( double ) * h5_time[bonds].size() );
        dataset_id = H5Dcreate(group_id, "time", H5T_NATIVE_DOUBLE, dataspace_id, 
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          h5feeder);
        status = H5Dclose(dataset_id);
        /*
        dataset_id = H5Dcreate(group_id, "px_last", H5T_NATIVE_DOUBLE, dataspace_id, 
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          h5_px_last[bonds]);
        status = H5Dclose(dataset_id);

        dataset_id = H5Dcreate(group_id, "px_bid", H5T_NATIVE_DOUBLE, dataspace_id, 
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          h5_px_bid[bonds]);
        status = H5Dclose(dataset_id);

        dataset_id = H5Dcreate(group_id, "px_low", H5T_NATIVE_DOUBLE, dataspace_id, 
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          h5_px_low[bonds]);
        status = H5Dclose(dataset_id);

        dataset_id = H5Dcreate(group_id, "px_high", H5T_NATIVE_DOUBLE, dataspace_id, 
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          h5_px_high[bonds]);
        status = H5Dclose(dataset_id);

        dataset_id = H5Dcreate(group_id, "px_ask", H5T_NATIVE_DOUBLE, dataspace_id, 
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          h5_px_ask[bonds]);
        status = H5Dclose(dataset_id);

        dataset_id = H5Dcreate(group_id, "bond_rate_last", H5T_NATIVE_DOUBLE, dataspace_id, 
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          h5_bond_rate_last[bonds]);
        status = H5Dclose(dataset_id);

        dataset_id = H5Dcreate(group_id, "bond_rate_bid", H5T_NATIVE_DOUBLE, dataspace_id, 
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          h5_bond_rate_bid[bonds]);
        status = H5Dclose(dataset_id);

        dataset_id = H5Dcreate(group_id, "bond_rate_low", H5T_NATIVE_DOUBLE, dataspace_id, 
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          h5_bond_rate_low[bonds]);
        status = H5Dclose(dataset_id);

        dataset_id = H5Dcreate(group_id, "bond_rate_high", H5T_NATIVE_DOUBLE, dataspace_id, 
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          h5_bond_rate_high[bonds]);
        status = H5Dclose(dataset_id);

        dataset_id = H5Dcreate(group_id, "bond_rate_ask", H5T_NATIVE_DOUBLE, dataspace_id, 
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          h5_bond_rate_ask[bonds]);
        status = H5Dclose(dataset_id);

        dataset_id = H5Dcreate(group_id, "zero_rate", H5T_NATIVE_DOUBLE, dataspace_id, 
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          zero_rate[bonds]);
        status = H5Dclose(dataset_id);
        */
        // Close the group
        status = H5Gclose(group_id);

        // Close the data space
        status = H5Sclose(dataspace_id);

        // Deallocate the native array (hdf5 feeder)
        delete[] h5feeder;

    };

    // Close the attribute data space
    status = H5Sclose(dataspace_id2);

    // Close the file
    status = H5Fclose(file_id);
};
