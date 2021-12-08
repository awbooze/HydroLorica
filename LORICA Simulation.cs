using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using LORICAVariables;


using System.Collections;
using System.Collections.Concurrent;
using System.Threading;
using System.Threading.Tasks;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.IO;
using System.Net;
using System.Runtime.InteropServices;
using System.Diagnostics;
using System.Numerics;
using MathNet.Numerics;
using MathNet.Numerics.IntegralTransforms;
using System.Windows.Forms;
using LORICA4;

namespace LORICA4
{
    class Simulation
    {
        Stopwatch stopwatch;
        GUIVariables guiVariables;

        int numfile;
        string filename;
        int maxruns, scenario, run_number;
        int save_time2, save_interval2 = 0;
        int coordinateDone = 0;
        double urfinalLati, urfinalLongi, llfinalLati, llfinalLongi, yurcorner, xurcorner = 0;
        TimeSpan geo_t, pedo_t, hydro_t, ponding_t;

        #region weird public variables
        AviWriter aw; // <JMW 20041018>
        bool frost_weathering_active,
            tilting_active,
            uplift_active,
            soil_phys_weath_active,
            soil_chem_weath_active,
            soil_bioturb_active,
            soil_clay_transloc_active,
            soil_carbon_active,
            memory_records,
            memory_records_d, 
            crashed,
            creep_active,
            water_ero_active,
            tillage_active,
            landslide_active,
            bedrock_weathering_active;

        int
                test,
                depressions_delta,
                depressions_alone,
                depressions_filled,  //counters for logging and reporting # of depressions filled/sedimented into/left alone
                depressionnumber = 0,
                maxdepressionnumber,
                depressionready,
                iloradius, iupradius,
                jloradius, jupradius,
                nbismemberofdepression,
                otherdepression,
                otherdepressionsize,
                totaldepressions,
                totaldepressionsize,
                maxsize,
                lowestneighbourcounter,
                numberoflowestneighbours,
                depressiondrainsout,
                largestdepression,
                xrow, xcol, xxrow, xxcol,
                landuse_value,
                number_of_outputs = 0,
                wet_cells, eroded_cells,
                i,
                j,
                low,
                high,
                equal,
                alpha,
                beta,
                numsinks,
                direct,
                deposited_cells;

        double soildepth_error, dtm00;
        //int GlobalMethods.n_texture_classes = 5;

        //int GlobalMethods.n_texture_classes = 5;

        const int maximumdepressionsize = 1500;  // run the program once to find out the SIZE OF the largest depression. Any higher number will do....
        const double tangent_of_delta = 0.005;
        const double tangent_for_outlet = 0.1;
        const int maxlowestnbs = 100000;
        const double epsilon = 0.000001;

        double[] local_s_i_t_kg = new double[] { 0, 0, 0, 0, 0 };

        // for constant layer thicknesses
        double dz_standard = 0.1;
        double tolerance = 0.55;

        int[] rainfall_record, evap_record, infil_record, till_record, temp_record;
        int[] rainfall_record_d, evap_record_d, duration_record_d;
        int[] zonesize = new int[22], zoneprogress = new int[22];

        int[]
        iloedge = new int[GlobalMethods.numberofsinks],
        jloedge = new int[GlobalMethods.numberofsinks],
        iupedge = new int[GlobalMethods.numberofsinks],
        jupedge = new int[GlobalMethods.numberofsinks],
        depressionsize = new int[GlobalMethods.numberofsinks],
        depressionconsidered = new int[GlobalMethods.numberofsinks],
        rowlowestnb = new int[maxlowestnbs],
        collowestnb = new int[maxlowestnbs];
        double available_for_delta_kg = 0;
        double available_for_delta_m = 0;

        //calibration globals
        int best_run, calib_levels, user_specified_number_of_calibration_parameters, user_specified_number_of_ratios;
        double reduction_factor, best_error;
        //USER INPUT NEEDED: establish best versions of parameters varied in calibration:
        double[] best_parameters;
        double[,] calib_ratios;
        double[] original_ratios;

        //soil timeseries_variables
        double total_average_soilthickness_m,
            total_phys_weathered_mass_kg,
            total_chem_weathered_mass_kg,
            total_fine_neoformed_mass_kg,
            total_fine_eluviated_mass_kg,
            total_mass_bioturbed_kg,
            total_OM_input_kg,
        local_soil_depth_m,
        local_soil_mass_kg;
        int number_soil_thicker_than,
        number_soil_coarser_than;



        // Water erosion and deposition parameters
        double
        advection_erodibility,
        P_act,
        m, n,				        // capacity slope and discharge exponents
        erosion_threshold_kg,
        rock_protection_constant,
        bio_protection_constant,
        constant_selective_transcap,
        Slope,			            // Gradient
        conv_fac,		            // convergence/divergence factor
        dS, desired_change, dztot,	// Difference in sediment/deposition/erosion
        sedtr_loc,                  // Local sediment transport rate
        all_grids,
        fraction,	                // fraction slope by slopesum
        frac_dis,	                // fraction of discharge into lower grid
        sediment_transported,		// fraction of transport rate
        water_out,
        unfulfilled_change, dz_left1, actual_change, 	// unfulfilled sedimentation
        dz_min, mmin,
        dz_max, maxx,		        // maximum lowest neighbour, steepest descent
        dz_bal, dz_bal2,		    // dz balans counter
        sedbal, sedbal2,
        erobal, erobal2,
        erobalto, sedbalto,
        erocnt,
        sedcnt,
        sediment_exported,		        // sediment out of our system
        total_Bolsena_sed_influx,

        // Biological weathering parameters  see Minasny and McBratney 2006 Geoderma 133
        P0,                         // m GlobalMethods.t-1  // weathering rate constant
        k1,                         // GlobalMethods.t-1
        k2,                         // GlobalMethods.t-1
        Pa,                         // m GlobalMethods.t-1  // weathering rate when soildepth = 0

        // Tilting and Uplift parameters
        tilt_intensity, lift_intensity,

        // Tillage parameter
        tilc;

        // Tree fall parameters
        double W_m_max, D_m_max, W_m, D_m, tf_frequency;
        int growth_a_max, age_a_max;

        //Soil physical weathering parameters
        double physical_weathering_constant, weathered_mass_kg, Cone, Ctwo;
        double[] upper_particle_size = new double[5];

        //Soil chemical weathering parameters
        double chemical_weathering_constant, Cthree, Cfour, Cfive, Csix, neoform_constant;
        double[] specific_area = new double[5];

        //Clay translocation parameters
        double max_eluviation, Cclay, ct_depthdec;

        //Bioturbation parameters
        double potential_bioturbation_kg;
        double bioturbation_depth_decay_constant;

        //Carbon cycle parameters
        double potential_OM_input,
               OM_input_depth_decay_constant,
               humification_fraction,
               potential_young_decomp_rate,
               potential_old_decomp_rate,
               young_depth_decay_constant,
               old_CTI_decay_constant,
               old_depth_decay_constant,
               young_CTI_decay_constant;

        // Decalcification parameters
        double[,,] CO3_kg;   // CaCO3, to track decalcification speed. Does not contribute to texture or soil mass (yet) MM
        double ini_CO3_content_frac;



        double noval,
        sediment_filled_m, depressionvolume_filled_m, sediment_delta_m,   // counters for logging and reporting the filling of depressions
        altidiff, minaltidiff,
        totaldepressionvolume,
        infil_value_m, evap_value_m, rain_value_m, soildepth_value,
        volume_eroded, volume_deposited,
        sum_normalweathered, sum_frostweathered, sum_soildepth, sum_creep, sum_solif, avg_solif, avg_creep, avg_soildepth,
        sum_ls, total_sum_tillage, total_sum_uplift, total_sum_tilting, total_sed_export;  // counters for logging and reporting through time

        int temp_value_C;

        double depressionsum_sediment_m, depressionsum_water_m, depressionsum_YOM_kg, depressionsum_OOM_kg;
        double[] depressionsum_texture_kg;
        double needed_to_fill_depression_m, dhoblique, dhobliquemax1, dhobliquemax2, firstalt, secondalt, dtmlowestnb;
        int dhmax_errors, readynum = 0, memberdepressionnotconsidered, depressionnum = 0, currentdepression;
        int lower_nb_exists, breaker = 0, rowlowestobnb, collowestobnb, II = 0, JJ = 0;
        int startrow, startcol, iloradius2, iupradius2, jupradius2, jloradius2, deltasize;
        int readysearching, iloradius3, iupradius3, jupradius3, jloradius3, couldbesink, omikron, omega;
        int tempx, tempy, obnbchanged;
        double sed_delta_size1 = 0, sed_delta_size2 = 0, sed_delta_size3 = 0;

        double
                diffusivity_creep,
                plough_depth,
                annual_weathering,
                dh, diff, dh1, dh_maxi,
                scan_do, dcount, powered_slope_sum,
                dmax, dmin,
                max_allowed_erosion,			// maximum erosion down to neighbour
                maximum_allowed_deposition, dhtemp,
                CSIZE,
                transport_capacity_kg,			// Capacity 
                detachment_rate,
                settlement_rate,
                frac_sed,   // fraction of landslide deposition into lower grids
                frac_bud,
                startsed,
                strsed,     // sediment delivered to streams
                T_act,         // Transmissivity
                bulkd_act,     // Bulk Density
                intfr_act,     // Internal Friction
                C_act,         // Combined Cohesion
                erotot,      // total landslide erosion
                sedtot,     // total landslide deposition;
                a_ifr, a_coh, a_bd, a_T,  // parameters parent material 1
                b_coh, b_ifr, b_bd, b_T,  // parameters parent material 2
                c_coh, c_ifr, c_bd, c_T,  // parameters parent material 3
                d_coh, d_ifr, d_bd, d_T,  // parameters parent material 4
                e_coh, e_ifr, e_bd, e_T,  // parameters parent material 5
                slopelim,       // slope limit for landslide erosion                          FACTOR 1
                celfrac,        // fraction used in calculation of celdistance (0.4 default)  FACTOR 2
                streamca,       // contributing area, number of cells, for stream development FACTOR 3
                rainfall_intensity,      // threshold critical rainfall value for landslide scenario   FACTOR 4
                slide_tot,
                dh_tot,
                tra_di,
                set_di,
                dh_tol,
                dt = 1,				// time step
                actual_t,      // Time counter for loop
                out_t,
                total_altitude,
                total_average_altitude,
                total_rain, total_evap, total_infil, total_outflow,
                //WVG
                total_sed_export_up, total_sed_export_mid, total_sed_export_low,
                total_sed_prod_up, total_sed_prod_mid, total_sed_prod_low,
                total_sed_dep_up, total_sed_dep_mid, total_sed_dep_low;  // counters for logging and reporting through time

        // tectonics
        int lift_type, lift_location, tilt_location;
        long scan_lon, scan_cnt, NRO, NCO;



        double[,] correct_dtm,        // for calibration purposes
                    paleo_dtm,          // for calibration purposes
                    sink_sed,
                    olddem,
                    dtm_epsilon,
                    timeseries_matrix,
                    profile_dtm1,       //WVG profile timeseries matrices
                    profile_dtm2,
                    profile_dtm3,
                    profile_wat1,
                    profile_wat2,
                    profile_wat3,
                    lessivage_errors; // for calibration of lessivage


        #endregion


        public Simulation(GUIVariables gv)
        {
            guiVariables = gv;


            depressionsum_texture_kg = new double[GlobalMethods.n_texture_classes];
        }


        public void main_loop(object sender, System.EventArgs e)
        {
            /*
            this.Invoke(new MethodInvoker(() => {
                InfoStatusPanel.Text = "Entered main program";
            }
            ));
            */
            guiVariables.InfoStatusPanel = "Entered main program";

            stopwatch = Stopwatch.StartNew();
            try
            {
                #region Fast Setup
                //foreach (string dtmfilename in Directory.EnumerateFiles(this.dtm_input_filename_textbox.Text, "*.txt", SearchOption.TopDirectoryOnly)) //"*.asc"
                //{
                string dtmfilename = guiVariables.DTM_input_filename_textbox;
                Debug.WriteLine("Entered LORICA main code with " + dtmfilename);
                string[] separate = dtmfilename.Split('.');
                GlobalMethods.Workdir = separate[0];
                Debug.WriteLine("storing results in " + GlobalMethods.Workdir);
                System.IO.Directory.CreateDirectory(GlobalMethods.Workdir);
                GlobalMethods.input_data_error = false;

                try { guiVariables.End_time = int.Parse(guiVariables.Number_runs_textbox); }
                catch { GlobalMethods.input_data_error = true; MessageBox.Show("Invalid number of years"); }
                /*try { ntr = System.Convert.ToInt32(guiVariables.End_time); }     // WVG initialise ntr: number of rows in timeseries matrix   
                catch (OverflowException)
                {
                    MessageBox.Show("number of timesteps is outside the range of the Int32 type.");
                }*/                                                                                                 //unused
                //WVG initialise ntr, GlobalMethods.nr of timesteps, can be changed to GlobalMethods.nr of output timesteps
                numfile = 1;
                guiVariables.ProcessStatusPanel = "";
                if (guiVariables.Water_ero_checkbox)
                {
                    water_ero_active = true; //Replace with guiVariables...
                    guiVariables.ProcessStatusPanel += "WE ";
                }
                if (guiVariables.Tillage_checkbox)
                {
                    tillage_active = true;  //Replace with guiVariables...
                    guiVariables.ProcessStatusPanel += "TI ";
                }
                if (guiVariables.Landslide_checkbox)
                {
                    landslide_active = true;  //Replace with guiVariables...
                    guiVariables.ProcessStatusPanel += "LS ";
                }
                if (guiVariables.Creep_active_checkbox)
                {
                    creep_active = true;  //Replace with guiVariables...
                    guiVariables.ProcessStatusPanel += "CR ";
                }
                if (guiVariables.Biological_weathering_checkbox)
                {
                    bedrock_weathering_active = true;  //Replace with guiVariables...
                    guiVariables.ProcessStatusPanel += "BW ";
                }
                if (guiVariables.Frost_weathering_checkbox)
                {
                    frost_weathering_active = true;  //Replace with guiVariables...
                    guiVariables.ProcessStatusPanel += "FW ";
                }
                if (guiVariables.Tilting_active_checkbox)
                {
                    tilting_active = true;  //Replace with guiVariables...
                    guiVariables.ProcessStatusPanel += "TL ";
                }
                if (guiVariables.Uplift_active_checkbox)
                {
                    uplift_active = true;  //Replace with guiVariables...
                    guiVariables.ProcessStatusPanel += "UP ";
                }
                if (guiVariables.Soil_phys_weath_checkbox)
                {
                    soil_phys_weath_active = true;  //Replace with guiVariables...
                    guiVariables.ProcessStatusPanel += "PW ";
                }
                if (guiVariables.Soil_chem_weath_checkbox)
                {
                    soil_chem_weath_active = true;  //Replace with guiVariables...
                    guiVariables.ProcessStatusPanel += "CW ";
                }
                if (guiVariables.Soil_bioturb_checkbox)
                {
                    soil_bioturb_active = true;  //Replace with guiVariables...
                    guiVariables.ProcessStatusPanel += "BT ";
                }
                if (guiVariables.Soil_clay_transloc_checkbox)
                {
                    soil_clay_transloc_active = true;  //Replace with guiVariables...
                    guiVariables.ProcessStatusPanel += "CT ";
                }
                if (guiVariables.Soil_carbon_cycle_checkbox) //:)
                {
                    soil_carbon_active = true;  //Replace with guiVariables...
                    guiVariables.ProcessStatusPanel += "CC ";
                }


                //INPUTS
                //GENERAL INPUTS
                //Entry point for consecutive runs for sensitivity analyses or calibration 
                maxruns = 1;
                int currentlevel = 0;

                if (guiVariables.Calibration_button)
                {
                    int runs_per_level = 0;
                    //CALIB_USER INPUT NEEDED NEXT LINE IN THE CODE :
                    user_specified_number_of_calibration_parameters = 1;
                    best_error = 99999999999; //or any other absurdly high number
                    best_parameters = new double[user_specified_number_of_calibration_parameters];
                    user_specified_number_of_ratios = guiVariables.Calibration_ratios_textbox.Split(';').Length;   //String made multiple times... Maybe single variable?
                    runs_per_level = Convert.ToInt32(Math.Pow(user_specified_number_of_ratios, user_specified_number_of_calibration_parameters));
                    calib_ratios = new double[user_specified_number_of_calibration_parameters, user_specified_number_of_ratios];
                    original_ratios = new double[user_specified_number_of_ratios];
                    for (int rat = 0; rat < user_specified_number_of_ratios; rat++)
                    {
                        try
                        {
                            original_ratios[rat] = Convert.ToDouble(guiVariables.Calibration_ratios_textbox.Split(';')[rat]);
                            for (int par = 0; par < user_specified_number_of_calibration_parameters; par++)
                            {
                                calib_ratios[par, rat] = Convert.ToDouble(guiVariables.Calibration_ratios_textbox.Split(';')[rat]);
                            }
                        }
                        catch { Debug.WriteLine(" problem setting original parameter ratios for calibration "); }
                    }
                    calib_calculate_maxruns(user_specified_number_of_calibration_parameters);
                    Debug.WriteLine(maxruns);
                    calib_prepare_report();
                    //CALIB_USER: set the number of parameters and their initial value

                }
                if (guiVariables.Sensitivity_button == true)
                { //dev needed
                }
                #endregion

                //parallelization possible here?
                for (run_number = 0; run_number < maxruns; run_number++)
                {
                    Debug.WriteLine(" maxruns is " + maxruns);

                    try { save_interval2 = System.Convert.ToInt32(guiVariables.GoogAnimationSaveInterval); }
                    catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for google save interval is not valid"); }

                    if (guiVariables.UTMgridcheckbox)
                    {
                        try { test = System.Convert.ToInt32(guiVariables.UTMzonebox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for UTM zone is not valid"); }
                    }

                    if (guiVariables.End_time < save_interval2 && guiVariables.GoogleAnimationCheckbox) { GlobalMethods.input_data_error = true; MessageBox.Show("value for google save interval cannot be larger than number of runs "); }
                    if (guiVariables.End_time < int.Parse(guiVariables.Saveintervalbox) && guiVariables.CheckBoxGenerateAVIFile) { GlobalMethods.input_data_error = true; MessageBox.Show("value for video interval cannot be larger than number of runs "); }


                    //WATER EROSION AND DEPOSITION PARAMETERS
                    if (water_ero_active)
                    {
                        try { m = double.Parse(guiVariables.Parameter_m_textbox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter m is not valid"); }                      // Kirkby's m and n factors for increasing
                        try { n = double.Parse(guiVariables.Parameter_n_textbox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter n is not valid"); }                   // sheet, wash, overland, gully to river flow
                        try { conv_fac = double.Parse(guiVariables.Parameter_m_textbox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter p is not valid"); }
                        try { advection_erodibility = double.Parse(guiVariables.Parameter_K_textbox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter K is not valid"); }
                        try { bio_protection_constant = double.Parse(guiVariables.Bio_protection_constant_textbox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter P is not valid"); }
                        try { rock_protection_constant = double.Parse(guiVariables.Rock_protection_constant_textbox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter P is not valid"); }
                        try { constant_selective_transcap = double.Parse(guiVariables.Selectivity_constant_textbox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter P is not valid"); }
                        try { erosion_threshold_kg = double.Parse(guiVariables.Erosion_threshold_textbox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter P is not valid"); }
                    }

                    //TILLAGE PARAMETERS
                    if (tillage_active)
                    {
                        try { plough_depth = double.Parse(guiVariables.Parameter_ploughing_depth_textbox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter plough depth is not valid"); }
                        try { tilc = double.Parse(guiVariables.Parameter_tillage_constant_textbox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter tillage constant is not valid"); }
                    }

                    //CREEP PARAMETER
                    if (creep_active)
                    {
                        try { conv_fac = double.Parse(guiVariables.Parameter_m_textbox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter p is not valid"); }
                        try { diffusivity_creep = double.Parse(guiVariables.Parameter_diffusivity_textbox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter diffusivity is not valid"); }
                    }

                    //LANDSLIDE PARAMETERS
                    if (landslide_active)
                    {
                        conv_fac = 4;        // multiple flow conversion factor
                    }

                    //Bio Weathering PARAMETERS
                    if (bedrock_weathering_active)
                    {
                        try { P0 = double.Parse(guiVariables.Parameter_P0_textbox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter P0 is not valid"); }
                        try { k1 = double.Parse(guiVariables.Parameter_k1_textbox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter k1 is not valid"); }
                        try { k2 = double.Parse(guiVariables.Parameter_k2_textbox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter k2 is not valid"); }
                        try { Pa = double.Parse(guiVariables.Parameter_Pa_textbox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter Pa is not valid"); }
                    }

                    //Tilting parameters
                    if (tilting_active)
                    {
                        if (guiVariables.Radio_tilt_col_zero) { tilt_location = 0; }
                        if (guiVariables.Radio_tilt_row_zero) { tilt_location = 1; }
                        if (guiVariables.Radio_tilt_col_max) { tilt_location = 2; }
                        if (guiVariables.Radio_tilt_row_max) { tilt_location = 3; }
                        try { tilt_intensity = double.Parse(guiVariables.Tilting_rate_textbox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter tilting rate is not valid"); }
                    }

                    //Uplift parameters
                    if (uplift_active)
                    {
                        if (guiVariables.Radio_lift_row_less_than) { lift_type = 0; }
                        if (guiVariables.Radio_lift_row_more_than) { lift_type = 1; }
                        if (guiVariables.Radio_lift_col_less_than) { lift_type = 2; }
                        if (guiVariables.Radio_lift_col_more_than) { lift_type = 3; }
                        if (lift_type == 0)
                        {
                            try { lift_location = int.Parse(guiVariables.Text_lift_row_less); }
                            catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter tilting rate is not valid"); }
                        }
                        if (lift_type == 1)
                        {
                            try { lift_location = int.Parse(guiVariables.Text_lift_row_more); }
                            catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter tilting rate is not valid"); }
                        }
                        if (lift_type == 2)
                        {
                            try { lift_location = int.Parse(guiVariables.Text_lift_col_less); }
                            catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter tilting rate is not valid"); }
                        }
                        if (lift_type == 3)
                        {
                            try { lift_location = int.Parse(guiVariables.Text_lift_col_more); }
                            catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter tilting rate is not valid"); }
                        }
                        try { lift_intensity = double.Parse(guiVariables.Uplift_rate_textbox); }
                        catch { GlobalMethods.input_data_error = true; MessageBox.Show("value for parameter tilting rate is not valid"); }
                    }

                    // TREE FALL PARAMETERS
                    if (guiVariables.Treefall_checkbox)
                    {
                        W_m_max = System.Convert.ToDouble(guiVariables.Tf_W);
                        D_m_max = System.Convert.ToDouble(guiVariables.Tf_D);
                        growth_a_max = System.Convert.ToInt32(guiVariables.Tf_growth);
                        age_a_max = System.Convert.ToInt32(guiVariables.Tf_age);
                        tf_frequency = System.Convert.ToDouble(guiVariables.Tf_freq);
                    }

                    //SOIL PHYSICAL WEATHERING PARAMETERS
                    if (soil_phys_weath_active)
                    {
                        try
                        {
                            physical_weathering_constant = Convert.ToDouble(guiVariables.Physical_weath_C1_textbox);
                            Cone = Convert.ToDouble(guiVariables.Physical_weath_constant1);
                            Ctwo = Convert.ToDouble(guiVariables.Physical_weath_constant2);
                            //the upper sizes of particle for the different fractions are declared in initialise_soil because they are always needed
                            Debug.WriteLine("succesfully read parameters for pysical weathering");
                        }
                        catch
                        {
                            GlobalMethods.input_data_error = true; Debug.WriteLine("problem reading parameters for pysical weathering");
                        }
                    }


                    //SOIL CHEMICAL WEATHERING PARAMETERS
                    if (soil_chem_weath_active)
                    {
                        try
                        {
                            chemical_weathering_constant = Convert.ToDouble(guiVariables.Chem_weath_rate_constant_textbox);
                            Cthree = Convert.ToDouble(guiVariables.Chem_weath_depth_constant_textbox);
                            Cfour = Convert.ToDouble(guiVariables.Chem_weath_specific_coefficient_textbox);
                            specific_area[0] = Convert.ToDouble(guiVariables.Specific_area_coarse_textbox);
                            specific_area[1] = Convert.ToDouble(guiVariables.Specific_area_sand_textbox);
                            specific_area[2] = Convert.ToDouble(guiVariables.Specific_area_silt_textbox);
                            specific_area[3] = Convert.ToDouble(guiVariables.Specific_area_clay_textbox);
                            specific_area[4] = Convert.ToDouble(guiVariables.Specific_area_fine_clay_textbox);
                            neoform_constant = Convert.ToDouble(guiVariables.Clay_neoform_constant_textbox);
                            Cfive = Convert.ToDouble(guiVariables.Clay_neoform_C1_textbox);
                            Csix = Convert.ToDouble(guiVariables.Clay_neoform_C2_textbox);
                            Debug.WriteLine("succesfully read parameters for chemical weathering");
                        }
                        catch
                        {
                            GlobalMethods.input_data_error = true; Debug.WriteLine("problem reading parameters for chemical weathering");
                        }
                    }

                    //SOIL CLAY DYNAMICS PARAMETERS
                    if (soil_clay_transloc_active)
                    {
                        try
                        {
                            max_eluviation = Convert.ToDouble(guiVariables.Maximum_eluviation_textbox);
                            Cclay = Convert.ToDouble(guiVariables.Eluviation_coefficient_textbox);
                            Debug.WriteLine("succesfully read parameters for  clay dynamics");
                        }
                        catch
                        {
                            GlobalMethods.input_data_error = true; Debug.WriteLine("problem reading parameters for clay dynamics");
                        }
                        if (guiVariables.CT_depth_decay_checkbox)
                        {
                            try
                            {
                                ct_depthdec = Convert.ToDouble(guiVariables.Ct_depth_decay);
                            }
                            catch
                            {
                                GlobalMethods.input_data_error = true; Debug.WriteLine("problem reading depth decay parameter for clay dynamics");
                            }
                        }
                    }

                    //BIOTURBATION PARAMETERS
                    if (soil_bioturb_active)
                    {
                        try
                        {
                            potential_bioturbation_kg = Convert.ToDouble(guiVariables.Potential_bioturbation_textbox);
                            bioturbation_depth_decay_constant = Convert.ToDouble(guiVariables.Bioturbation_depth_decay_textbox);
                        }
                        catch
                        {
                            GlobalMethods.input_data_error = true; Debug.WriteLine("problem reading parameters for bioturbation");
                        }
                    }

                    //CARBON CYCLE PARAMETERS
                    if (soil_carbon_active)
                    {
                        try
                        {
                            potential_OM_input = Convert.ToDouble(guiVariables.Carbon_input_textbox);
                            OM_input_depth_decay_constant = Convert.ToDouble(guiVariables.Carbon_depth_decay_textbox);
                            humification_fraction = Convert.ToDouble(guiVariables.Carbon_humification_fraction_textbox);
                            potential_young_decomp_rate = Convert.ToDouble(guiVariables.Carbon_y_decomp_rate_textbox);
                            potential_old_decomp_rate = Convert.ToDouble(guiVariables.Carbon_o_decomp_rate_textbox);
                            young_depth_decay_constant = Convert.ToDouble(guiVariables.Carbon_y_depth_decay_textbox);
                            old_CTI_decay_constant = Convert.ToDouble(guiVariables.Carbon_o_twi_decay_textbox);
                            old_depth_decay_constant = Convert.ToDouble(guiVariables.Carbon_o_depth_decay_textbox);
                            young_CTI_decay_constant = Convert.ToDouble(guiVariables.Carbon_y_twi_decay_textbox);
                        }
                        catch
                        {
                            GlobalMethods.input_data_error = true; Debug.WriteLine("problem reading parameters for carbon cycle");
                        }
                    }

                    try
                    {
                        filename = dtmfilename;             //for directory input
                        GlobalMethods.dtm_file(filename);                 // from dtm_file(), almost all memory for the model is claimed
                    }
                    catch { Debug.WriteLine(" failed to initialise GlobalMethods.dtm "); }

                    //LARGEST THING IN HERE. RUNS FOREVER (very long time)
                    if (!GlobalMethods.input_data_error)
                    {
                        try
                        {

                            //Debug.WriteLine("reading general values");
                            if (!guiVariables.Check_space_soildepth)
                            {
                                try { soildepth_value = double.Parse(guiVariables.Soildepth_constant_value_box); }
                                catch { MessageBox.Show("value for parameter soildepth is not valid"); }
                            }
                            if (!guiVariables.Check_space_landuse && !guiVariables.Check_time_landuse)
                            {
                                try { landuse_value = int.Parse(guiVariables.Landuse_constant_value_box); }
                                catch { MessageBox.Show("value for parameter GlobalMethods.landuse is not valid"); }
                            }
                            if (!guiVariables.Check_space_evap && !guiVariables.Check_time_evap)
                            {
                                try { evap_value_m = double.Parse(guiVariables.Evap_constant_value_box); }
                                catch { MessageBox.Show("value for parameter GlobalMethods.evapotranspiration is not valid"); }
                            }
                            if (!guiVariables.Check_space_infil && !guiVariables.Check_time_infil)
                            {
                                try { infil_value_m = double.Parse(guiVariables.Infil_constant_value_box); }
                                catch { MessageBox.Show("value for parameter infiltration is not valid"); }
                            }

                            if (!guiVariables.Check_space_rain && !guiVariables.Check_time_rain)
                            {
                                try { rain_value_m = double.Parse(guiVariables.Rainfall_constant_value_box); }
                                catch { MessageBox.Show("value for parameter rainfall is not valid"); }
                            }

                            if (!guiVariables.Check_time_T)
                            {
                                try { temp_value_C = int.Parse(guiVariables.Temp_constant_value_box); }
                                catch { MessageBox.Show("value for parameter temperature is not valid"); }
                            }

                        }
                        catch { MessageBox.Show("there was a problem reading input values"); GlobalMethods.input_data_error = true; }
                        // Debug.WriteLine("initialising non-general inputs");
                        try { initialise_once(); } // reading input files
                        catch { MessageBox.Show("there was a problem reading input files "); GlobalMethods.input_data_error = true; }


                        //CALIB_USER: multiply parameter values with current ratio
                        //Note the correspondence between the formulas. Change only 1 value for additional parameters!
                        if (guiVariables.Calibration_button)
                        {
                            //Debug.WriteLine("erodib " + advection_erodibility + " conv fac " + conv_fac);
                            int rat_number = Convert.ToInt32(Math.Floor(run_number / Math.Pow(user_specified_number_of_ratios, 0)) % user_specified_number_of_ratios);
                            advection_erodibility *= calib_ratios[0, rat_number];
                            Debug.WriteLine("First ratio number: " + rat_number);
                            //rat_number = Convert.ToInt32(Math.Floor(run_number / Math.Pow(user_specified_number_of_ratios, 1)) % user_specified_number_of_ratios);
                            //conv_fac *= calib_ratios[1, rat_number];
                            //Debug.WriteLine("Second ratio number: " + rat_number);
                            // Debug.WriteLine("erodib " + advection_erodibility + " conv fac " + conv_fac);
                        }

                        timeseries_matrix = new double[System.Convert.ToInt32(guiVariables.End_time), number_of_outputs];
                        if (GlobalMethods.input_data_error == false)
                        {
                            try
                            {
                                /*tabControl1.Visible = false;  
                                Mapselector.Enabled = true; 
                                try {View_tabs_checkbox.Checked = false;}
                                catch { Debug.WriteLine(" failed to set view_tabs to unchecked "); }
                                graphicToGoogleEarthButton.Visible = true;
                                graphics_scale = 4;
                                double c_scale = System.Convert.ToDouble(780) / System.Convert.ToDouble(GlobalMethods.nc);
                                double r_scale = System.Convert.ToDouble(330) / System.Convert.ToDouble(GlobalMethods.nr);
                                if (c_scale < r_scale) { graphics_scale = System.Convert.ToInt32(Math.Floor(c_scale)); }
                                else { graphics_scale = System.Convert.ToInt32(Math.Floor(r_scale)); }
                                if (graphics_scale < 1) { graphics_scale = 1; }
                                m_objDrawingSurface = new Bitmap(GlobalMethods.nc * graphics_scale, GlobalMethods.nr * graphics_scale, System.Drawing.Imaging.PixelFormat.Format24bppRgb);
                                mygraphics = this.CreateGraphics();
                                Mapwindow.Visible = true;
                                //map_controls.Visible = true; */
                            }
                            catch { Debug.WriteLine("graphics initialisation failed "); GlobalMethods.input_data_error = true; }
                            if (GlobalMethods.input_data_error == false)
                            {
                                int count_intervene = 0;
                                GlobalMethods.t_intervene = 0;
                                if (GlobalMethods.t_intervene > 0) { read_soil_elevation_distance_from_output(GlobalMethods.t_intervene, GlobalMethods.Workdir); }

                                //begining of looping
                                for (GlobalMethods.t = GlobalMethods.t_intervene; GlobalMethods.t < guiVariables.End_time; GlobalMethods.t++)
                                {

                                    try
                                    {
                                        every_timestep();
                                    }
                                    catch
                                    {
                                        Debug.WriteLine("failed to run in timestep " + GlobalMethods.t);
                                        // Catch for when the model crashes due to unknown reasons. The model will read the latest output and start calculating again from there which I named an intervention). When the crash occurs five times, the model breaks MM
                                        if (count_intervene < 5)
                                        {
                                            count_intervene += 1;
                                            GlobalMethods.t_intervene = GlobalMethods.t - (GlobalMethods.t % (int.Parse(guiVariables.Box_years_output)));
                                            Debug.WriteLine("intervening at GlobalMethods.t" + GlobalMethods.t_intervene);
                                            read_soil_elevation_distance_from_output(GlobalMethods.t_intervene, GlobalMethods.Workdir);
                                        }
                                        else
                                        {
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        MessageBox.Show("input data error - program can not yet run");
                        //tabControl1.Visible = true;
                        //Mapselector.Enabled = false;
                        //map_controls.Visible = false;
                    }

                    if (guiVariables.Calibration_button)
                    {
                        //calculate how good this run was:
                        double current_error = calib_objective_function();
                        //store that information along with the parameter values used to achieve it:
                        calib_update_report(current_error);
                        if (current_error < best_error) { best_error = current_error; calib_update_best_paras(); best_run = run_number; }
                        //and check whether one 'level' of calibration has finished. If so, we have to change parameter values
                        Debug.WriteLine("run " + run_number + " number paras " + user_specified_number_of_calibration_parameters + " number ratios " + guiVariables.Calibration_ratios_textbox.Split(';').Length);
                        if ((run_number + 1) % Convert.ToInt32(Math.Pow(guiVariables.Calibration_ratios_textbox.Split(';').Length, user_specified_number_of_calibration_parameters)) == 0)
                        {

                            //a level of calibration has finished

                            //If it was the last level, we are now done
                            currentlevel++;
                            Debug.WriteLine(" successfully finished a level of calibration runs");
                            if (run_number == maxruns - 1)
                            {
                                Debug.WriteLine(" successfully finished last level of calibration runs");
                                calib_finish_report();
                            }
                            else
                            {
                                Debug.WriteLine(" setting new ratios ");
                                //CALIB_USER INPUT NEEDED HERE IN CODE
                                //check whether the best run was on the edge of parameter space or inside, shift to that place and zoom out or in
                                calib_shift_and_zoom(0, double.Parse(guiVariables.Calibration_ratio_reduction_parameter_textbox), double.Parse(guiVariables.Parameter_K_textbox));
                                //calib_shift_and_zoom(1, double.Parse(calibration_ratio_reduction_parameter_textbox.Text), double.Parse(guiVariables.Parameter_m_textbox));
                            }
                        }
                        else
                        {
                            //nothing. Parameter values are adapted with the corresponding ratios to continue calibration above.
                        }
                    }

                } // exit point for consecutive runs


            } // end try
            catch
            {
                Debug.WriteLine("Error in accessing file " + guiVariables.DTM_input_filename_textbox);
            }

            guiVariables.UpdateAllFields();

        }  //end main

        //start of running code


        private void calculate_overwater_landscape_Spitsbergen()
        {
            //to account for a landscape that is isostaically rebounding from below sealevel to above sealevel. 
            //height above sealevel itself is not important, just that the landscape grows over time
            //therefore, Marijn's solution: if (elevation < threshold(GlobalMethods.t)) , then elevation = nodata
            double minimum_overwater_elevation = (10263 - GlobalMethods.t) / 218;
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    if (GlobalMethods.original_dtm[GlobalMethods.row, GlobalMethods.col] != -9999)
                    {

                        if (GlobalMethods.original_dtm[GlobalMethods.row, GlobalMethods.col] > minimum_overwater_elevation && GlobalMethods.dtm[GlobalMethods.row, GlobalMethods.col] == -9999)
                        {
                            // these cases were not over water, but now will be.
                            GlobalMethods.dtm[row, col] = GlobalMethods.original_dtm[row, col];
                            // all cases that were already overwater, will stay overwater - no changes there.
                        }
                    }
                    else
                    {
                        // nothing happens because these cells are simply not part of the study area
                    }
                }
            }
        }

        private void every_timestep()    //performs actions in every timestep
        {
            // If a cell should remain at the same fixed elevation (e.g. fixed elevation boundary condition), here the cell can be selected
            //if(GlobalMethods.nr>50&GlobalMethods.nc>100)
            //{ 
            //    if (GlobalMethods.t == GlobalMethods.t_intervene) { dtm00 = GlobalMethods.dtm[50, 100]; }

            //    // force no-change boundary op de outlet CLORPT
            //    GlobalMethods.dtm[50, 100] = dtm00;

            //}

            DateTime geo_start, pedo_start, hydro_start;


            //if (GlobalMethods.t == 0 | GlobalMethods.t == 1) { displaysoil(50, 0); }
            int i = 0;
            this.TimeStatusPanel.Text = "timestep " + (GlobalMethods.t + 1) + "/" + +guiVariables.End_time;
            updateClick = 1;
            // Debug.WriteLine("starting calculations - TIME " + GlobalMethods.t);

            if (GlobalMethods.Ik_ben_Marijn)
            { calculate_overwater_landscape_Spitsbergen(); }
            //displaysoil(0, 0);

            #region hydrological processes
            //Debug.WriteLine("before water {0}",GlobalMethods.texture_kg[0,0,0,2]);
            //Debug.WriteLine("before water balance");
            hydro_start = DateTime.Now;
            if (guiVariables.Daily_water)
            {

                // Debug.WriteLine("Running daily water balance");
                water_balance();
                // Debug.WriteLine("Daily water balance finished");

            }
            //print_spatial_water_balance();
            // print_P_ET0();
            hydro_t += DateTime.Now - hydro_start;
            #endregion

            #region Vegetation
            // Debug.WriteLine("before vegetation");
            if (guiVariables.Daily_water)
            {
                determine_vegetation_type();
                change_vegetation_parameters();
            }

            #endregion

            #region Geomorphic processes
            geo_start = DateTime.Now;

            //Debug.WriteLine("before WE");
            //displaysoil(0, 0);
            if (water_ero_active)
            {
                //Debug.WriteLine("before WE1");

                initialise_every(); //fast
                GlobalMethods.comb_sort();        //fast

                if (guiVariables.Daily_water)
                {
                    //Debug.WriteLine("before WE2");
                    calculate_water_ero_sed_daily();
                    //Debug.WriteLine("before WE3");
                    soil_update_split_and_combine_layers();
                    //Debug.WriteLine("before WE4");

                }
                else
                {
                    //Debug.WriteLine("calculating water erosion");
                    findsinks();
                    searchdepressions();
                    define_fillheight_new();
                    if (NA_anywhere_in_soil() == true) { Debug.WriteLine("NA found before erosed"); }
                    calculate_water_ero_sed();
                    //Debug.WriteLine("Overland flow in GlobalMethods.t " + GlobalMethods.t+" with a flow of "+Math.Round(rain_value_m-infil_value_m-evap_value_m,4)+" m");
                    if (NA_anywhere_in_soil() == true) { Debug.WriteLine("NA found after erosed"); }
                    if (crashed) { Debug.WriteLine("crashed while calculating water erosion"); }
                    //else { Debug.WriteLine("successfully finished water erosion calculation"); }

                }
            }

            // Debug.WriteLine("before TF");
            if (guiVariables.Treefall_checkbox)
            {
                if (GlobalMethods.t <= (guiVariables.End_time - 500)) // if there is no tillage
                {
                    calculate_tree_fall();
                }

                // Debug.WriteLine("Calculating tree fall");
            }

            //displaysoil(0, 0);

            if (bedrock_weathering_active)
            {
                //Debug.WriteLine("calculating bedrock weathering");
                calculate_bedrock_weathering();
                soil_update_split_and_combine_layers();
                //if (GlobalMethods.t % 25 == 0) { displaysoil(55, 108); }
            }
            // Debug.WriteLine("before CR");
            //displaysoil(0, 0);
            if (creep_active)
            {
                try
                {// Debug.WriteLine("calculating GlobalMethods.creep");
                    GlobalMethods.comb_sort();

                    calculate_creep();

                    soil_update_split_and_combine_layers();
                }
                catch { Debug.WriteLine(" failed during GlobalMethods.creep calculations"); }
            }


            // Debug.WriteLine("before TI");
            //displaysoil(0, 0);
            if (tillage_active)
            {
                GlobalMethods.comb_sort();
                int tilltime = 0;
                //if (guivariables.Check_time_till_fields) { tilltime = till_record[GlobalMethods.t]; }
                //else { tilltime = 1; }

                if (GlobalMethods.t > (guiVariables.End_time - 500)) { tilltime = 1; }
                if (tilltime == 1)
                {


                    initialise_every_till();
                    calculate_tillage();
                    soil_update_split_and_combine_layers();
                }

            }
            //Debug.WriteLine("after TI");
            //displaysoil(0, 0);
            if (landslide_active)
            {
                Debug.WriteLine("calculating landsliding");
                GlobalMethods.comb_sort();
                ini_slope();
                calculate_critical_rain();
                calculate_slide();
            }

            geo_t += DateTime.Now - geo_start;

            #endregion

            #region pedogenic processes


            pedo_start = DateTime.Now;

            // Debug.WriteLine("before PW");
            //displaysoil(0, 0);
            if (soil_phys_weath_active)
            {
                // Debug.WriteLine("calculating soil physical weathering");
                if (!GlobalMethods.Ik_ben_Marijn) { soil_physical_weathering(); }
                else
                {
                    SPITS_soil_physical_weathering();
                    SPITS_aeolian_deposition();
                }
                soil_update_split_and_combine_layers();


            }
            // Debug.WriteLine("before CW");
            //displaysoil(0, 0);
            if (soil_chem_weath_active)
            {
                //Debug.WriteLine("calculating soil chemical weathering");
                soil_chemical_weathering();
                soil_update_split_and_combine_layers();
                if (guiVariables.Timeseries.Total_average_soilthickness_checkbox)
                {
                    timeseries_matrix[GlobalMethods.t, timeseries_order[21]] = total_average_soilthickness_m;
                }
                if (guiVariables.Timeseries.Timeseries_number_soil_thicker_checkbox)
                {
                    timeseries_matrix[GlobalMethods.t, timeseries_order[22]] = number_soil_thicker_than;
                }
                if (guiVariables.Timeseries.Timeseries_coarser_checkbox)
                {
                    timeseries_matrix[GlobalMethods.t, timeseries_order[23]] = number_soil_coarser_than;
                }
                if (guiVariables.Timeseries.Timeseries_soil_depth_checkbox)
                {
                    timeseries_matrix[GlobalMethods.t, timeseries_order[24]] = local_soil_depth_m;
                }
                if (guiVariables.Timeseries.Timeseries_soil_mass_checkbox)
                {
                    timeseries_matrix[GlobalMethods.t, timeseries_order[25]] = local_soil_mass_kg;
                }
            }

            // Debug.WriteLine("before CT");
            //displaysoil(0, 0);
            if (soil_clay_transloc_active)
            {
                // Debug.WriteLine("calculating soil clay dynamics ");

                if (GlobalMethods.Ik_ben_Marijn)
                {
                    soil_silt_translocation(); // Spitsbergen case study
                }
                else
                {
                    if (guiVariables.Ct_Jagercikova)
                    {
                        soil_clay_translocation_Jagercikova();
                    }
                    else
                    {
                        soil_clay_translocation();
                    }
                }
                soil_update_split_and_combine_layers();
                if (NA_in_map(GlobalMethods.dtm) > 0 | NA_in_map(GlobalMethods.soildepth_m) > 0)
                {
                    Debug.WriteLine("err_ets1");
                }

            }
            if (NA_anywhere_in_soil() == true) { Debug.WriteLine("NA found before soil carbon"); }
            //displaysoil(0, 0);
            if (soil_carbon_active)
            {
                // Debug.WriteLine("calculating carbon dynamics ");
                if (guiVariables.Version_lux_checkbox)
                {
                    soil_litter_cycle();
                }
                else
                {
                    soil_carbon_cycle();
                }

            }
            if (NA_anywhere_in_soil() == true) { Debug.WriteLine("NA found after soil carbon"); }
            if (guiVariables.Decalcification_checkbox)
            {
                //Debug.WriteLine("calculating decalcification");
                soil_decalcification();
            }

            if (soil_bioturb_active)
            {
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        update_all_soil_thicknesses(row, col);
                    }
                }
                // Debug.WriteLine("calculating bioturbation");
                soil_bioturbation();
                // if (findnegativetexture()) { Debugger.Break(); }

                soil_update_split_and_combine_layers();
                // if (findnegativetexture()) { Debugger.Break(); }

            }
            if (NA_anywhere_in_soil() == true) { Debug.WriteLine("NA found after soil bioturb"); }


            pedo_t += DateTime.Now - pedo_start;

            #endregion

            #region write output

            // Debug.WriteLine("before output");
            //if (view_maps_checkbox.Checked == true && GlobalMethods.t == guiVariables.End_time - 1) { draw_map(mygraphics); updateClick = 1; }
            if (guiVariables.View_maps_checkbox == true && GlobalMethods.t == guiVariables.End_time - 1) { updateClick = 1; }
            numfile++;





            int t_out = GlobalMethods.t + 1;
            if ((guiVariables.Final_output_checkbox && t_out == guiVariables.End_time) || (guiVariables.Regular_output_checkbox && ((t_out) % (int.Parse(guiVariables.Box_years_output)) == 0)))
            {
                if (GlobalMethods.t == guiVariables.End_time - 1)
                {

                    //Debug.WriteLine("Time balance. Geomorphic processes: {0} min, pedogenic processes: {1} min, hydrologic processes: {2} min, ponding {3} min", geo_t, pedo_t, hydro_t, ponding_t);
                }
                //Debug.WriteLine("Attempting to write outputs");

                // displaysoil(31, 12);
                // Debug.WriteLine("Total catchment mass = " + total_catchment_mass());

                if (guiVariables.Daily_water)
                {
                    Debug.WriteLine("writing daily water");

                    //try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_aridity.asc", aridity_vegetation); }
                    //catch { MessageBox.Show("vegetation has not been written"); }

                    try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_infiltration_m.asc", Iy); }
                    catch { MessageBox.Show("infiltration has not been written"); }

                    try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_actual_evapotranspiration_m.asc", ETay); }
                    catch { MessageBox.Show("ETa has not been written"); }

                    try
                    {
                        GlobalMethods.out_integer(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_vegetationtype.asc", vegetation_type);
                        
                        for (int row = 0; row < GlobalMethods.nr; row++)
                        {
                            for (int col = 0; col < GlobalMethods.nc; col++)
                            {
                                GlobalMethods.vegetation_type[row, col] = 0; // reset vegetation_type, to give the output per output period
                            }
                        }
                    }
                    catch { MessageBox.Show("vegetation type has not been written"); }



                }

                if (guiVariables.Version_lux_checkbox)
                {
                    try
                    {
                        // outputs for case study Luxembourg. Focus on different litter types
                        // young labile OM is hornbeam, old stable OM is beech
                        // Outputs:
                        // SOM stocks entire profile: total, young, old (kg/m2)
                        // top layer: total, old, young (-) 

                        string[] litter_types = { "hornbeam", "beech", "total" };
                        string[] litter_outputs = { "stocks_kgm2", "toplayer_frac" };

                        foreach (string type in litter_types) // loop over different SOM types
                        {
                            // determine which SOM fraction should be considered
                            bool h_bool = false; bool b_bool = false;
                            if (type == "hornbeam") { h_bool = true; }
                            if (type == "beech") { b_bool = true; }
                            if (type == "total") { h_bool = true; b_bool = true; }

                            foreach (string output in litter_outputs)
                            {
                                // determine which layers to consider and what to calculate
                                int numberoflayers = 0;
                                if (output == "stocks_kgm2") { numberoflayers = GlobalMethods.max_soil_layers; }
                                if (output == "toplayer_frac") { numberoflayers = 1; }

                                double[,] output_litter_map = new double[GlobalMethods.nr, GlobalMethods.nc];

                                for (int row = 0; row < GlobalMethods.nr; row++)
                                {
                                    for (int col = 0; col < GlobalMethods.nc; col++)
                                    {
                                        double litterstock_kg = 0;
                                        double mineralsoil_toplayer_kg = 0;
                                        if (h_bool) { litterstock_kg += GlobalMethods.litter_kg[row, col, 0]; }
                                        if (b_bool) { litterstock_kg += GlobalMethods.litter_kg[row, col, 1]; }

                                        if (output == "toplayer_frac")
                                        {
                                            for (int tex = 0; tex < 5; tex++)
                                            {
                                                mineralsoil_toplayer_kg += GlobalMethods.texture_kg[row, col, 0, tex];
                                            }
                                        }


                                        if (output == "toplayer_frac") { litterstock_kg /= (mineralsoil_toplayer_kg + litterstock_kg); } // calculate to fraction
                                        if (output == "stocks_kgm2") { litterstock_kg /= (GlobalMethods.dx * GlobalMethods.dx); } // calculate to kg/m2
                                        output_litter_map[row, col] = litterstock_kg;
                                    }
                                }
                                try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_litter_" + type + "_" + output + ".asc", output_litter_map); }
                                catch { MessageBox.Show("litter output has not been written"); }
                            }
                        }
                        try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_TPI.asc", GlobalMethods.tpi); }
                        catch { MessageBox.Show("TPI output has not been written"); }

                        /* CODE BLOCK BELOW WRITES OUT DIFFERENT ORGANIC MATTER MAPS. THIS IS NOT NECESSARY ANYMORE NOW LITTER IS STORED IN ITS OWN MATRIX
                         * 
                        // outputs for case study Luxembourg. Focus on different organic matter types
                        // young labile OM is hornbeam, old stable OM is beech
                        // Outputs:
                        // SOM stocks entire profile: total, young, old (kg/m2)
                        // top layer: total, old, young (-) 

                        string[] SOM_types = { "young", "old", "total" };
                        string[] SOM_outputs = { "stocks_kgm2", "toplayer_frac" };

                        foreach (string type in SOM_types) // loop over different SOM types
                        {
                            // determine which SOM fraction should be considered
                            bool y_bool = false; bool o_bool = false;
                            if (type == "young") { y_bool = true; }
                            if (type == "old") { o_bool = true; }
                            if (type == "total") { y_bool = true; o_bool = true; }

                            foreach (string output in SOM_outputs)
                            {
                                // determine which layers to consider and what to calculate
                                int numberoflayers = 0;
                                if (output == "stocks_kgm2") { numberoflayers = GlobalMethods.max_soil_layers; }
                                if (output == "toplayer_frac") { numberoflayers = 1; }

                                double[,] output_SOM_map = new double[GlobalMethods.nr, GlobalMethods.nc] ;
                                for (int row = 0; row < GlobalMethods.nr; row++)
                                {
                                    for (int col = 0; col < GlobalMethods.nc; col++)
                                    {
                                        double SOMstock_kg = 0;
                                        double mineralsoil_kg = 0;
                                        for (int lay = 0; lay < numberoflayers; lay++)
                                        {
                                            if (y_bool) { SOMstock_kg += GlobalMethods.young_SOM_kg[row, col, lay]; }
                                            if (o_bool) { SOMstock_kg += GlobalMethods.old_SOM_kg[row, col, lay]; }

                                            if (output == "toplayer_frac")
                                            {
                                                for (int tex = 0; tex < 5; tex++)
                                                {
                                                    mineralsoil_kg += GlobalMethods.texture_kg[row, col, lay, tex];
                                                }
                                            }
                                        }
                                        if (output == "toplayer_frac") { SOMstock_kg /= (mineralsoil_kg + SOMstock_kg); } // calculate to fraction
                                        if (output == "stocks_kgm2") { SOMstock_kg /= (GlobalMethods.dx * GlobalMethods.dx); } // calculate to kg/m2
                                        output_SOM_map[row, col] = SOMstock_kg;
                                    }
                                }
                                try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_SOM_" + type + "_" + output + ".asc", output_SOM_map); }
                                catch { MessageBox.Show("SOM output has not been written"); }
                            }
                        }
                        try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_TPI.asc", GlobalMethods.tpi); }
                        catch { MessageBox.Show("TPI output has not been written"); }

                        */
                    }
                    catch
                    {
                        Debug.WriteLine("Error in writing litterwater_ outputs for Luxembourg case study");
                    }
                }



                try
                {
                    //Debug.WriteLine("writing all soils");
                    GlobalMethods.writeallsoils();
                }
                catch
                {
                    Debug.WriteLine("Failed during writing of soils");
                }


                if (guiVariables.Altitude_output_checkbox)
                {

                    try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_dtm.asc", GlobalMethods.dtm); }
                    catch { MessageBox.Show("dtm has not been written"); }

                    try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_dz_soil.asc", GlobalMethods.dz_soil); }
                    catch { MessageBox.Show("dz_soil has not been written"); }


                    //try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + GlobalMethods.t + "_out_dzero.asc", GlobalMethods.dz_ero_m); }
                    //catch { MessageBox.Show("dzero has not been written"); }
                    //try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + GlobalMethods.t + "_out_dzsed.asc", GlobalMethods.dz_sed_m); }
                    //catch { MessageBox.Show("dzsed has not been written"); }
                }
                if (guiVariables.Treefall_checkbox)
                {
                    try
                    {
                        GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_dz_treefall.asc", GlobalMethods.dz_treefall);
                        GlobalMethods.out_integer(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_treefallcount.asc", GlobalMethods.treefall_count);

                    }
                    catch { MessageBox.Show("treefall has not been written"); }
                }
                if (guiVariables.Soildepth_output_checkbox)
                {
                    try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_soildepth.asc", GlobalMethods.soildepth_m); }
                    catch { MessageBox.Show("soildepth has not been written"); }
                }
                if (guiVariables.Alt_change_output_checkbox)
                {
                    try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_change.asc", GlobalMethods.dtmchange); }
                    catch { MessageBox.Show("change has not been written"); }
                }

                if (guiVariables.Water_output_checkbox & guiVariables.Water_ero_checkbox)
                {
                    // Debug.WriteLine("before writing water flow");


                    try
                    {
                        if (guiVariables.Daily_water)
                        {
                            for (int roww = 0; roww < GlobalMethods.nr; roww++)
                            {
                                for (int colw = 0; colw < GlobalMethods.nc; colw++)
                                {
                                    GlobalMethods.waterflow_m3[roww, colw] = guiVariables.OFy_m[roww, colw, 0];
                                }
                            }
                        }
                        GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_water.asc", GlobalMethods.waterflow_m3);
                    }
                    catch { MessageBox.Show("water has not been written"); }
                }
                if (guiVariables.Depressions_output_checkbox)
                {
                    try { GlobalMethods.out_integer(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_depress.asc", depression); }
                    catch { MessageBox.Show("depressions have not been written"); }
                    try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_dtmfillA.asc", GlobalMethods.dtmfill_A); }
                    catch { MessageBox.Show("dfmfill has not been written"); }
                }
                if (guiVariables.Diagnostic_output_checkbox)
                {
                    //try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + GlobalMethods.t + "_out_sedintrans.asc", sediment_in_transport); }
                    //catch {  MessageBox.Show("sed in trans has not been written"); }
                    try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_dzero.asc", GlobalMethods.dz_ero_m); }
                    catch { MessageBox.Show("dzero has not been written"); }
                    try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_dzsed.asc", GlobalMethods.dz_sed_m); }
                    catch { MessageBox.Show("dzsed has not been written"); }
                    try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_lakesed.asc", GlobalMethods.lake_sed_m); }
                    catch { MessageBox.Show("lakesed has not been written"); }
                }

                if (guiVariables.Water_ero_checkbox)
                {
                    // Debug.WriteLine("before writing water erosion");

                    if (guiVariables.All_process_output_checkbox)
                    {
                        try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_water_erosion.asc", GlobalMethods.sum_water_erosion); }
                        catch { MessageBox.Show("water erosion has not been written"); }
                    }
                }
                if (creep_active_checkbox.Checked)
                {
                    // Debug.WriteLine("before writing GlobalMethods.creep");

                    try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_creep.asc", GlobalMethods.creep); }
                    catch { MessageBox.Show("GlobalMethods.creep has not been written"); }

                }

                if (Tillage_checkbox.Checked)
                {
                    // Debug.WriteLine("before writing tillage erosion");

                    if (all_process_output_checkbox.Checked)
                    {
                        try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_out_tillage.asc", GlobalMethods.sum_tillage); }
                        catch { MessageBox.Show("tillage has not been written"); }
                    }
                }
                if (guiVariables.Landslide_checkbox)
                {
                    try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_GlobalMethods.crrain.asc", GlobalMethods.crrain); }
                    catch { MessageBox.Show("GlobalMethods.crrain has not been written"); }
                    try { GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_ca.asc", GlobalMethods.camf); }
                    catch { MessageBox.Show("ca has not been written"); }
                }

                if (guiVariables.Decalcification_checkbox)
                {
                    try
                    {
                        double[,] decalcification_depth = new double[GlobalMethods.nr, GlobalMethods.nc];
                        for (int rowdec = 0; rowdec < GlobalMethods.nr; rowdec++)
                        {
                            for (int coldec = 0; coldec < GlobalMethods.nc; coldec++)
                            {
                                bool decal_written = false;
                                if (GlobalMethods.dtm[rowdec, coldec] != -9999)
                                {
                                    double depthdec = 0;
                                    int laydec = 0;
                                    while (decal_written == false)
                                    {
                                        if (CO3_kg[rowdec, coldec, laydec] == 0 && (laydec != (GlobalMethods.max_soil_layers - 1)))
                                        {
                                            if (laydec < (GlobalMethods.max_soil_layers - 1))
                                            {
                                                laydec++;
                                                depthdec += GlobalMethods.layerthickness_m[rowdec, coldec, laydec];
                                            }

                                        }
                                        else
                                        {
                                            decalcification_depth[rowdec, coldec] = depthdec;
                                            decal_written = true;
                                        }
                                    }
                                }
                            }
                        }
                        GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + t_out + "_decaldepth.asc", decalcification_depth);
                    }
                    catch
                    {
                        MessageBox.Show("decalcification has not been written");
                    }
                }

                if (guiVariables.Profile.Radio_pro1_col)
                {
                    if (guiVariables.Profile.Check_altitude_profile1)
                    {
                        try { GlobalMethods.out_profile(GlobalMethods.Workdir + "\\profile_1_dtm_" + run_number + "_" + t_out + ".asc", GlobalMethods.dtm, false, System.Convert.ToInt32(guiVariables.Profile.P1_row_col_box)); }
                        catch { MessageBox.Show("profile_1_dtm_" + run_number + "_" + t_out + ".asc has not been written"); }
                    }
                    if (guiVariables.Profile.Check_waterflow_profile1)
                    {
                        try { GlobalMethods.out_profile(GlobalMethods.Workdir + "\\profile_1_water_" + run_number + "_" + t_out + ".asc", GlobalMethods.waterflow_m3, false, System.Convert.ToInt32(guiVariables.Profile.P1_row_col_box)); }
                        catch { MessageBox.Show("profile_1_water_" + run_number + "_" + t_out + ".asc has not been written"); }
                    }
                }
                if (guiVariables.Profile.Radio_pro1_row)
                {
                    if (guiVariables.Profile.Check_altitude_profile1)
                    {
                        try { GlobalMethods.out_profile(GlobalMethods.Workdir + "\\profile_1_dtm_" + run_number + "_" + t_out + ".asc", GlobalMethods.dtm, true, System.Convert.ToInt32(guiVariables.Profile.P1_row_col_box)); }
                        catch { MessageBox.Show("profile_1_dtm_" + run_number + "_" + t_out + ".asc has not been written"); }
                    }
                    if (guiVariables.Profile.Check_waterflow_profile1)
                    {
                        try { GlobalMethods.out_profile(GlobalMethods.Workdir + "\\profile_1_water_" + run_number + "_" + t_out + ".asc", GlobalMethods.waterflow_m3, true, System.Convert.ToInt32(guiVariables.Profile.P1_row_col_box)); }
                        catch { MessageBox.Show("profile_1_water_" + run_number + "_" + t_out + ".asc has not been written"); }
                    }
                }
                if (guiVariables.Profile.Radio_pro2_col)
                {
                    if (guiVariables.Profile.Check_altitude_profile1)
                    {
                        try { GlobalMethods.out_profile(GlobalMethods.Workdir + "\\profile_2_dtm_" + run_number + "_" + t_out + ".asc", GlobalMethods.dtm, false, System.Convert.ToInt32(guiVariables.Profile.P2_row_col_box)); }
                        catch { MessageBox.Show("profile_2_dtm_" + run_number + "_" + t_out + ".asc has not been written"); }
                    }
                    if (guiVariables.Profile.Check_waterflow_profile1)
                    {
                        try { GlobalMethods.out_profile(GlobalMethods.Workdir + "\\profile_2_water_" + run_number + "_" + t_out + ".asc", GlobalMethods.waterflow_m3, false, System.Convert.ToInt32(guiVariables.Profile.P2_row_col_box)); }
                        catch { MessageBox.Show("profile_2_water_" + run_number + "_" + t_out + ".asc has not been written"); }
                    }
                }
                if (guiVariables.Profile.Radio_pro2_row)
                {
                    if (guiVariables.Profile.Check_altitude_profile1)
                    {
                        try { GlobalMethods.out_profile(GlobalMethods.Workdir + "\\profile_2_dtm_" + run_number + "_" + t_out + ".asc", GlobalMethods.dtm, true, System.Convert.ToInt32(guiVariables.Profile.P2_row_col_box)); }
                        catch { MessageBox.Show("profile_dtm_" + run_number + "_" + t_out + ".asc has not been written"); }
                    }
                    if (guiVariables.Profile.Check_waterflow_profile1)
                    {
                        try { GlobalMethods.out_profile(GlobalMethods.Workdir + "\\profile_2_water_" + run_number + "_" + t_out + ".asc", GlobalMethods.waterflow_m3, true, System.Convert.ToInt32(guiVariables.Profile.P2_row_col_box)); }
                        catch { MessageBox.Show("profile_2_water_" + run_number + "_" + t_out + ".asc has not been written"); }
                    }
                }
                if (guiVariables.Profile.Radio_pro3_col)
                {
                    if (guiVariables.Profile.Check_altitude_profile1)
                    {
                        try { GlobalMethods.out_profile(GlobalMethods.Workdir + "\\profile_3_dtm_" + run_number + "_" + t_out + ".asc", GlobalMethods.dtm, false, System.Convert.ToInt32(guiVariables.Profile.P3_row_col_box)); }
                        catch { MessageBox.Show("profile_3_dtm_" + run_number + "_" + t_out + ".asc has not been written"); }
                    }
                    if (guiVariables.Profile.Check_waterflow_profile1)
                    {
                        try { GlobalMethods.out_profile(GlobalMethods.Workdir + "\\profile_3_water_" + run_number + "_" + t_out + ".asc", GlobalMethods.waterflow_m3, false, System.Convert.ToInt32(guiVariables.Profile.P3_row_col_box)); }
                        catch { MessageBox.Show("profile_3_water_" + run_number + "_" + t_out + ".asc has not been written"); }
                    }
                }
                if (guiVariables.Profile.Radio_pro3_row)
                {
                    if (guiVariables.Profile.Check_altitude_profile1)
                    {
                        try { GlobalMethods.out_profile(GlobalMethods.Workdir + "\\profile_3_dtm_" + run_number + "_" + t_out + ".asc", GlobalMethods.dtm, true, System.Convert.ToInt32(guiVariables.Profile.P3_row_col_box)); }
                        catch { MessageBox.Show("profile_3_dtm_" + run_number + "_" + t_out + ".asc has not been written"); }
                    }
                    if (guiVariables.Profile.Check_waterflow_profile1)
                    {
                        try { GlobalMethods.out_profile(GlobalMethods.Workdir + "\\profile_3_water_" + run_number + "_" + t_out + ".asc", GlobalMethods.waterflow_m3, true, System.Convert.ToInt32(guiVariables.Profile.P3_row_col_box)); }
                        catch { MessageBox.Show("profile_3_water_" + run_number + "_" + t_out + ".asc has not been written"); }
                    }
                }
                //Debug.WriteLine("after outputs");

            }
            //Google Earth Animation
            if ((googleAnimationCheckbox.Checked) && (GlobalMethods.t % (int.Parse(googAnimationSaveInterval.Text)) == 0))
            {
                try { Google_Earth_Output(); }
                catch { Debug.WriteLine("Error writing Google Output"); }
            }
            //AVI file output
            if ((checkBoxGenerateAVIFile.Checked) && (GlobalMethods.t % (int.Parse(saveintervalbox.Text)) == 0))
            {
                try { AVI_Output(); }
                catch { Debug.WriteLine("Error writing video Output"); }
            }
            //if (GlobalMethods.t % 10 == 1)
            //{
            //    //displaysoil(25, 50);
            //}
            if (GlobalMethods.t == guiVariables.End_time - 1)
            {


                try
                {
                    //close google earth animation
                    if (googleAnimationCheckbox.Checked == true)
                    {
                        StreamWriter kmlsr = File.AppendText(KML_FILE_NAME);
                        kml = "\n</Folder>"
                              + "\n</kml>";
                        kmlsr.WriteLine(kml);
                        kmlsr.Close();
                    }
                }
                catch
                {
                    Debug.WriteLine("Error finishing Google Earth Output");
                }
                //close AVI file
                try
                {
                    if (checkBoxGenerateAVIFile.Checked)
                        aw.Close();  //JMW 20041109
                }
                catch (Exception ex)
                {
                    Debug.WriteLine("Error finishing video output");
                }

                this.InfoStatusPanel.Text = " --finished--";
                stopwatch.Stop();
                Debug.WriteLine("Elapsed time: " + stopwatch.Elapsed);
                //Timeseries output
                if (number_of_outputs > 0) { timeseries_output(); }
            }
            #endregion

        }

        private void Google_Earth_Output()
        {
            updateClick = 1;
            this.Refresh();
            draw_map(mygraphics);

            if (coordinateDone == 0)
            {
                //transfrom coordinates
                point testPoint = new point(GlobalMethods.xcoord, GlobalMethods.ycoord);
                if (UTMgridcheckbox.Checked)
                {
                    testPoint.UTMzone = System.Convert.ToInt32(UTMzonebox.Text);
                    testPoint.south = System.Convert.ToBoolean(UTMsouthcheck.Checked);
                    testPoint.transformUTMPoint();
                }
                else
                {
                    testPoint.transformPoint();
                }
                yurcorner = GlobalMethods.ycoord + (System.Convert.ToDouble(GlobalMethods.nr) * System.Convert.ToDouble(GlobalMethods.dx)); //ART possibly incorrect GlobalMethods.nr = GlobalMethods.nc
                xurcorner = GlobalMethods.xcoord + (System.Convert.ToDouble(GlobalMethods.nc) * System.Convert.ToDouble(GlobalMethods.dx));
                point testPoint2 = new point(xurcorner, yurcorner);
                if (UTMgridcheckbox.Checked)
                {
                    testPoint2.UTMzone = System.Convert.ToInt32(UTMzonebox.Text);
                    testPoint2.south = System.Convert.ToBoolean(UTMsouthcheck.Checked);
                    testPoint2.transformUTMPoint();
                }
                else
                {
                    testPoint2.transformPoint();
                }
                urfinalLati = testPoint2.ycoord;
                urfinalLongi = testPoint2.xcoord;
                llfinalLati = testPoint.ycoord;
                llfinalLongi = testPoint.xcoord;
                coordinateDone = 1;
            }

            //Save image
            m_objDrawingSurface.MakeTransparent();
            m_objDrawingSurface.Save(GlobalMethods.Workdir + "\\animation\\mysavedimage" + imageCount2 + ".png", System.Drawing.Imaging.ImageFormat.Png);
            //update time
            googleTime = googleTime.AddYears(save_interval2);
            kmlTime = googleTime.ToString();
            DateArray = kmlTime.Split(new char[] { ' ' });
            DateArray2 = DateArray[0].Split(new char[] { '-' });
            kmlTime = DateArray2[2] + "-" + DateArray2[1] + "-" + DateArray2[0] + "T" + DateArray[1] + "Z";

            //create kml file for image
            StreamWriter kmlsr = File.AppendText(KML_FILE_NAME);
            if (imageCount2 == 1)
            {
                kml = @"<?xml version=""1.0"" encoding=""UTF-8""?>
                         <kml xmlns=""http://earth.google.com/kml/2.1"">";
                kml = kml + "\n<Folder>"
                    + "\n<name>Animation</name>";
                kmlsr.WriteLine(kml);
                kml = "";
            }
            kml = kml + "\n<GroundOverlay>"
                + "\n<name>Untitled Image Overlay</name>";
            kml = kml + "\n<TimeSpan>"
                   + "\n<begin>" + kmlTime + "</begin>"
                   + "\n<end>" + kmlTime + "</end>"
                   + "\n</TimeSpan>"
                   + "\n<Icon>"
                   + "\n<href>mySavedImage" + imageCount2 + ".png</href>"
                   + "\n</Icon>"
                   + "\n<LatLonBox>";
            kml = kml + "\n<north>" + urfinalLati + "</north>"
                  + "\n<south>" + llfinalLati + "</south>"
                  + "\n<east>" + urfinalLongi + "</east>"
                  + "\n<west>" + llfinalLongi + "</west>\n";
            kml = kml + @"</LatLonBox>
                           </GroundOverlay>";
            kmlsr.WriteLine(kml);
            kml = "";
            kmlsr.Close();
            imageCount2 = imageCount2 + 1;

        }

        private void AVI_Output()
        {
            this.Refresh(); // tjc to enable graphics to be drawn before sending to AVI
            draw_map(mygraphics); // tjc
            Graphics gbmp = Graphics.FromImage(bmp);

            if (gbmp != null)
            {

                IntPtr dc1 = mygraphics.GetHdc();
                IntPtr dc2 = gbmp.GetHdc();

                //BitBlt(dc2, 0, 0, this.ClientRectangle.Width, this.ClientRectangle.Height,
                //dc1, 0, 0, 13369376);
                // this makes sure the entire LORICA window gets video-ed.

                BitBlt(dc2, this.Mapwindow.Location.X, this.Mapwindow.Location.Y,
                    this.Mapwindow.Size.Width + this.Mapwindow.Location.X,
                    this.Mapwindow.Size.Height + this.Mapwindow.Location.Y,
                    dc1, 0, 0, 13369376);
                // this makes sure the video only gets made for the mapped area


                mygraphics.ReleaseHdc(dc1);
                gbmp.ReleaseHdc(dc2);


                // need to flip image to get it the correct way up in the avi - not sure why.
                bmp.RotateFlip(System.Drawing.RotateFlipType.RotateNoneFlipY);


                try
                {
                    aw.AddFrame();
                }
                catch (AviWriter.AviException ex)  // <JMW 20041018>
                {
                    aw.Close();
                    Debug.WriteLine("AVI Exception in: " + ex.ToString());
                }
                try
                {
                    aw.AddFrame();
                }
                catch (AviWriter.AviException ex)  // <JMW 20041018>
                {
                    aw.Close();
                    Debug.WriteLine("AVI Exception in: " + ex.ToString());
                }
                try
                {
                    aw.AddFrame();
                }
                catch (AviWriter.AviException ex)  // <JMW 20041018>
                {
                    aw.Close();
                    Debug.WriteLine("AVI Exception in: " + ex.ToString());
                }
            }
        }

        private void read_soil_elevation_distance_from_output(int time, string dir)
        {
            // read latest output and start calculating from there
            dir = dir + "\\";

            initialise_once();

            filename = dir + "0_" + time + "_out_dtm.asc";
            GlobalMethods.read_double(filename, GlobalMethods.dtm);
            Debug.WriteLine("read GlobalMethods.dtm");

            filename = dir + "0_" + time + "_out_soildepth.asc";
            GlobalMethods.read_double(filename, GlobalMethods.soildepth_m);
            Debug.WriteLine("read soildepth");

            filename = dir + "0_" + time + "_out_change.asc";
            GlobalMethods.read_double(filename, GlobalMethods.dtmchange);
            Debug.WriteLine("read GlobalMethods.dtm change");

            filename = dir + "0_" + time + "_out_water_erosion.asc";
            GlobalMethods.read_double(filename, GlobalMethods.sum_water_erosion);
            Debug.WriteLine("read water erosion");

            filename = dir + "0_" + time + "_out_tillage.asc";
            GlobalMethods.read_double(filename, GlobalMethods.sum_tillage);
            Debug.WriteLine("read GlobalMethods.sum_tillage");

            filename = dir + "0_" + time + "_out_creep.asc";
            GlobalMethods.read_double(filename, GlobalMethods.creep);
            Debug.WriteLine("read GlobalMethods.creep");

            filename = dir + "0_" + time + "_out_dz_treefall.asc";
            GlobalMethods.read_double(filename, GlobalMethods.dz_treefall);
            Debug.WriteLine("read GlobalMethods.dz_treefall");

            filename = dir + "0_" + time + "_out_tillage.asc";
            GlobalMethods.read_double(filename, GlobalMethods.sum_tillage);
            Debug.WriteLine("read GlobalMethods.sum_tillage");

            filename = dir + "0_" + time + "_out_dz_soil.asc";
            GlobalMethods.read_double(filename, GlobalMethods.dz_soil);
            Debug.WriteLine("read sum_dz_soil");
            // SOIL INFORMATION
            // reset old info
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    for (int lay = 0; lay < GlobalMethods.max_soil_layers; lay++)
                    {
                        GlobalMethods.texture_kg[row, col, lay, 0] = 0;
                        GlobalMethods.texture_kg[row, col, lay, 1] = 0;
                        GlobalMethods.texture_kg[row, col, lay, 2] = 0;
                        GlobalMethods.texture_kg[row, col, lay, 3] = 0;
                        GlobalMethods.texture_kg[row, col, lay, 4] = 0;
                        GlobalMethods.young_SOM_kg[row, col, lay] = 0;
                        GlobalMethods.old_SOM_kg[row, col, lay] = 0;
                        GlobalMethods.layerthickness_m[row, col, lay] = 0;
                        GlobalMethods.bulkdensity[row, col, lay] = 0;
                    }
                }
            }


            using (var reader = new StreamReader(dir + "GlobalMethods.t" + time + "_out_allsoils.csv"))
            {
                int row, col, lay;
                // discard first line (header)    
                var line = reader.ReadLine();
                var values = line.Split(',');

                while (!reader.EndOfStream)
                {
                    line = reader.ReadLine();
                    values = line.Split(',');

                    row = Convert.ToInt32(values[0]);
                    col = Convert.ToInt32(values[1]);
                    lay = Convert.ToInt32(values[3]);

                    GlobalMethods.texture_kg[row, col, lay, 0] = Convert.ToDouble(values[8]); //coarse
                    GlobalMethods.texture_kg[row, col, lay, 1] = Convert.ToDouble(values[9]); // sand
                    GlobalMethods.texture_kg[row, col, lay, 2] = Convert.ToDouble(values[10]); // silt
                    GlobalMethods.texture_kg[row, col, lay, 3] = Convert.ToDouble(values[11]); // clay
                    GlobalMethods.texture_kg[row, col, lay, 4] = Convert.ToDouble(values[12]); // fine clay
                    GlobalMethods.young_SOM_kg[row, col, lay] = Convert.ToDouble(values[13]); // young SOM
                    GlobalMethods.old_SOM_kg[row, col, lay] = Convert.ToDouble(values[14]); // old SOM
                    GlobalMethods.layerthickness_m[row, col, lay] = Convert.ToDouble(values[5]); // thickness
                    GlobalMethods.bulkdensity[row, col, lay] = Convert.ToDouble(values[23]); // bulk density

                }
            }
        }


        private void timeseries_output()
        {
            int step;
            string FILENAME = workdir + "\\timeseries.log";
            using (StreamWriter sw = new StreamWriter(FILENAME))
            {
                //geomprph centred
                if (guiVariables.Timeseries.Timeseries_cell_waterflow_check) { sw.Write("cell_waterflow "); }
                if (guiVariables.Timeseries.Timeseries_cell_altitude_check) { sw.Write("cell_altitude "); }
                if (guiVariables.Timeseries.Timeseries_net_ero_check) { sw.Write("net_erosion "); }
                if (guiVariables.Timeseries.Timeseries_number_dep_check) { sw.Write("deposited_cells "); }
                if (guiVariables.Timeseries.Timeseries_number_erosion_check) { sw.Write("eroded_cells "); }
                if (guiVariables.Timeseries.Timeseries_number_waterflow_check) { sw.Write("wet_cells "); }
                if (guiVariables.Timeseries.Timeseries_SDR_check) { sw.Write("SDR "); }
                if (guiVariables.Timeseries.Timeseries_total_average_alt_check) { sw.Write("average_alt "); }
                if (guiVariables.Timeseries.Timeseries_total_dep_check) { sw.Write("total_dep "); }
                if (guiVariables.Timeseries.Timeseries_total_ero_check) { sw.Write("total_ero "); }
                if (guiVariables.Timeseries.Timeseries_total_evap_check) { sw.Write("total_evap "); }
                if (guiVariables.Timeseries.Timeseries_total_infil_check) { sw.Write("total_infil "); }
                if (guiVariables.Timeseries.Timeseries_total_outflow_check) { sw.Write("total_outflow "); }
                if (guiVariables.Timeseries.Timeseries_total_rain_check) { sw.Write("total_rain "); }
                //soil_centred
                if (guiVariables.Timeseries.Total_phys_weath_checkbox) { sw.Write("phys_weath "); }
                if (guiVariables.Timeseries.Total_chem_weath_checkbox) { sw.Write("chem_weath "); }
                if (guiVariables.Timeseries.Total_fine_formed_checkbox) { sw.Write("fine_clay_formed "); }
                if (guiVariables.Timeseries.Total_fine_eluviated_checkbox) { sw.Write("fine_clay_eluviated "); }
                if (guiVariables.Timeseries.Total_mass_bioturbed_checkbox) { sw.Write("mass_bioturbed "); }
                if (guiVariables.Timeseries.Total_OM_input_checkbox) { sw.Write("OM_input "); }
                if (guiVariables.Timeseries.Total_average_soilthickness_checkbox) { sw.Write("average_soilthickness "); }
                if (guiVariables.Timeseries.Timeseries_number_soil_thicker_checkbox) { sw.Write("soil_thicker "); }
                if (guiVariables.Timeseries.Timeseries_number_soil_thicker_checkbox) { sw.Write("soil_coarser "); }
                if (guiVariables.Timeseries.Timeseries_number_soil_thicker_checkbox) { sw.Write("soil_thickness "); }
                if (guiVariables.Timeseries.Timeseries_number_soil_thicker_checkbox) { sw.Write("soil_mass "); }
                sw.Write("\r\n");
                for (step = 0; step <= guiVariables.End_time - 1; step++)
                {
                    if (guiVariables.Timeseries.Timeseries_cell_waterflow_check) { sw.Write(timeseries_matrix[step, timeseries_order[1]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Timeseries_cell_altitude_check) { sw.Write(timeseries_matrix[step, timeseries_order[2]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Timeseries_net_ero_check) { sw.Write(timeseries_matrix[step, timeseries_order[3]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Timeseries_number_dep_check) { sw.Write(timeseries_matrix[step, timeseries_order[4]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Timeseries_number_erosion_check) { sw.Write(timeseries_matrix[step, timeseries_order[5]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Timeseries_number_waterflow_check) { sw.Write(timeseries_matrix[step, timeseries_order[6]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Timeseries_SDR_check) { sw.Write(timeseries_matrix[step, timeseries_order[7]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Timeseries_total_average_alt_check) { sw.Write(timeseries_matrix[step, timeseries_order[8]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Timeseries_total_dep_check) { sw.Write(timeseries_matrix[step, timeseries_order[9]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Timeseries_total_ero_check) { sw.Write(timeseries_matrix[step, timeseries_order[10]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Timeseries_total_evap_check) { sw.Write(timeseries_matrix[step, timeseries_order[11]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Timeseries_total_infil_check) { sw.Write(timeseries_matrix[step, timeseries_order[12]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Timeseries_total_outflow_check) { sw.Write(timeseries_matrix[step, timeseries_order[13]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Timeseries_total_rain_check) { sw.Write(timeseries_matrix[step, timeseries_order[14]]); sw.Write(" "); }
                    //soil_centred
                    if (guiVariables.Timeseries.Total_phys_weath_checkbox) { sw.Write(timeseries_matrix[step, timeseries_order[15]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Total_chem_weath_checkbox) { sw.Write(timeseries_matrix[step, timeseries_order[16]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Total_fine_formed_checkbox) { sw.Write(timeseries_matrix[step, timeseries_order[17]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Total_fine_eluviated_checkbox) { sw.Write(timeseries_matrix[step, timeseries_order[18]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Total_mass_bioturbed_checkbox) { sw.Write(timeseries_matrix[step, timeseries_order[19]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Total_OM_input_checkbox) { sw.Write(timeseries_matrix[step, timeseries_order[20]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Total_average_soilthickness_checkbox) { sw.Write(timeseries_matrix[step, timeseries_order[21]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Timeseries_number_soil_thicker_checkbox) { sw.Write(timeseries_matrix[step, timeseries_order[22]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Timeseries_coarser_checkbox) { sw.Write(timeseries_matrix[step, timeseries_order[23]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Timeseries_soil_depth_checkbox) { sw.Write(timeseries_matrix[step, timeseries_order[24]]); sw.Write(" "); }
                    if (guiVariables.Timeseries.Timeseries_soil_mass_checkbox) { sw.Write(timeseries_matrix[step, timeseries_order[25]]); }
                    sw.Write("\r\n");
                }
            }
        }


        #region Hydrology code

        //int[] guiVariables.P_all, ET0_all, guiVariables.Tavg_all, guiVariables.Tmin_all, guiVariables.Tmax_all, guiVariables.D_all;
        double[] Py = new double[365], ET0_m = new double[12];
        int[] Tavgy = new int[365], Tminy = new int[365], Tmaxy = new int[365], Dy = new int[365];
        // int[] guiVariables.D_all = new int[123], Dy = new int[365];
        double[,,] Ks_md, water_balance_m, Ra_rcm /*, guiVariables.OFy_m*/;
        double[,] Iy, ROy, Ks_topsoil_mh, pond_d, pond_y, outflow_y, stagdepth, waterfactor, total_outflow_y, ETay, ET0y;
        int[] month = new int[12] { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
        int[] monthcum = new int[12] { 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 };
        int[] midmonthdays = new int[] { 16, 46, 75, 106, 136, 167, 197, 228, 259, 289, 320, 350 };
        double Ks_min_mh, Ks_max_mh, snow_m, snow_start_m, snowfall_m, snowmelt_factor_mTd, snow_threshold_C, Pduringsnowmelt; // snow thickness is now constant in space. develop: spatially varying effects of snowfall. snowmelt factor in m T-1 d-1
                                                                                                                               // snowmelt modeling in Hock, 2003, Eq. 2: https://www.sciencedirect.com/science/article/pii/S0022169403002579#BIB50
                                                                                                                               // read data into guiVariables.P_all  etc
        DateTime ponding_start;

        double local_solar_radiation(double slope_rad, double aspect_rad, int month)
        {
            //Debug.WriteLine("lsr1");
            // SOURCE: Swift 1976: Algorithm for solar radiation on mountain slopes
            // https://doi.org/10.1029/WR012i001p00108
            double lat_rad, L1_Ra, L2_Ra, D1_Ra, D_Ra, E_Ra, R0_Ra, R1_Ra, R4_Ra, T_Ra, T7_Ra, T6_Ra, T3_Ra, T2_Ra, T1_Ra, T0_Ra, acos_in;

            R0_Ra = 1.95 * 0.041868; // Convert from cal/cm2/min to MJ/m2/min
            lat_rad = Math.PI / 180 * (System.Convert.ToDouble(latitude_deg.Text) + System.Convert.ToDouble(latitude_min.Text) / 60); // latitude [rad]
            L1_Ra = Math.Asin(Math.Cos(slope_rad) * Math.Sin(lat_rad) + Math.Sin(slope_rad) * Math.Cos(lat_rad) * Math.Cos(aspect_rad)); // equivalent latitude [rad]
            D1_Ra = Math.Cos(slope_rad) * Math.Cos(lat_rad) - Math.Sin(slope_rad) * Math.Sin(lat_rad) * Math.Cos(aspect_rad);
            if (D1_Ra == 0) { D1_Ra = 1E-10; }
            L2_Ra = Math.Atan(Math.Sin(slope_rad) * Math.Sin(aspect_rad) / D1_Ra);
            if (D1_Ra < 0) { L2_Ra += Math.PI; }

            int day_Ra = midmonthdays[month];

            D_Ra = 0.4 * Math.PI / 180 - 23.3 * Math.PI / 180 * Math.Cos((day_Ra + 10) * Math.PI / 180 * 0.986);
            E_Ra = 1 - 0.0167 * Math.Cos((day_Ra - 3) * Math.PI / 180 * 0.986);
            R1_Ra = 60 * R0_Ra / (E_Ra * E_Ra);

            acos_in = -Math.Tan(L1_Ra) * Math.Tan(D_Ra);
            if (acos_in > 1) { acos_in = 1; }
            T_Ra = Math.Acos(acos_in);
            T7_Ra = T_Ra - L2_Ra;
            T6_Ra = -T_Ra - L2_Ra;

            acos_in = -Math.Tan(lat_rad) * Math.Tan(D_Ra);
            if (acos_in > 1) { acos_in = 1; }
            T_Ra = Math.Acos(acos_in);
            T1_Ra = T_Ra;
            T0_Ra = -T_Ra;
            if (T7_Ra < T1_Ra) { T3_Ra = T7_Ra; } else { T3_Ra = T1_Ra; }
            if (T6_Ra > T0_Ra) { T2_Ra = T6_Ra; } else { T2_Ra = T0_Ra; }

            R4_Ra = R1_Ra * (Math.Sin(D_Ra) * Math.Sin(L1_Ra) * (T3_Ra - T2_Ra) * 12 / Math.PI + Math.Cos(D_Ra) * Math.Cos(L1_Ra) * (Math.Sin(T3_Ra + L2_Ra) - Math.Sin(T2_Ra + L2_Ra)) * 12 / Math.PI);
            //Debug.WriteLine("lsr2");
            if (R4_Ra < 0)
            {
                Debug.WriteLine("err_lsr1");
            }
            return (R4_Ra * 0.408 / 1000); // convert to m/d
        }

        void update_solar_radiation()
        {
            if (GlobalMethods.t % 100 == 0)
            {
                GlobalMethods.update_slope_and_aspect(); // updates slopemap GlobalMethods.slopeAnalysis [rad] and GlobalMethods.aspect [rad]
            }

            for (int hrow = 0; hrow < GlobalMethods.nr; hrow++)
            {
                for (int hcol = 0; hcol < GlobalMethods.nc; hcol++)
                {
                    if (GlobalMethods.dtm[hrow, hcol] != -9999)
                    {
                        for (int mo = 0; mo < 12; mo++)
                        {
                            //Debug.WriteLine("sr1");
                            // Fill 3D matrix with local monthly ET
                            Ra_rcm[hrow, hcol, mo] = local_solar_radiation(GlobalMethods.slopeAnalysis[hrow, hcol], GlobalMethods.aspect[hrow, hcol], mo);
                            //Debug.WriteLine("sr2");
                        }

                    }
                }
            }
        }

        double total_snow_melt, total_water_flow;
        void water_balance()
        {
            // Debug.WriteLine("wb1");
            total_water_flow = 0;
            if (GlobalMethods.t % 100 == 0)
            {
                update_solar_radiation();
            }
            Pduringsnowmelt = 0;

            total_snow_melt = 0;
            snow_start_m = snow_m;
            snowfall_m = 0;
            create_daily_weather(); // Calculate daily weather variables

            // Create yearly matrices for infiltration and overland flow
            // Debug.WriteLine("wb2");
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    pond_y[row, col] = 0;
                    outflow_y[row, col] = 0;
                    waterfactor[row, col] = 1;
                    total_outflow_y[row, col] = 0;
                    ETay[row, col] = 0;
                    ET0y[row, col] = 0;
                    Iy[row, col] = 0;
                    for (int dir = 0; dir < 10; dir++)
                    {
                        guiVariables.OFy_m[row, col, dir] = 0;
                    }
                }
            }
            GlobalMethods.comb_sort();

            double P_OF, snowmelt_m = 0;
            // Debug.WriteLine("wb3");
            if (GlobalMethods.t % 10 == 0) { update_Ks(); }

            // Debug.WriteLine("wb4");
            int daycount = 0;
            for (int mo = 0; mo < 12; mo++)
            {
                // create monthly timeseries
                double[] Pm = new double[month[mo]];
                int[] Dm = new int[month[mo]], Tavgm = new int[month[mo]];
                P_OF = 0;
                Array.Copy(Py, daycount, Pm, 0, month[mo]);
                Array.Copy(Dy, daycount, Dm, 0, month[mo]);
                Array.Copy(Tavgy, daycount, Tavgm, 0, month[mo]);
                ;
                // daily overland flow
                snow_threshold_C = 0;
                for (int day = 0; day < month[mo]; day++)
                {
                    if (Tavgm[day] <= snow_threshold_C) // temperature below 0, all GlobalMethods.rain falls as snow
                    {
                        if (Pm[day] > 0)
                        {
                            // Debug.WriteLine("snowfall");
                            snowfall_m += Pm[day];
                            snow_m += Pm[day];
                            P_OF += Pm[day];
                        }
                    }
                    else // T above 0, snow can melt and is all added to runoff. 
                    {
                        if (snow_m > 0) // Snow present
                        {
                            Pduringsnowmelt += Pm[day];
                            snowmelt_m = snowmelt_factor_mTd * (Tavgm[day] - snow_threshold_C);
                            if (snowmelt_m > snow_m) { snowmelt_m = snow_m; }
                            snow_m -= snowmelt_m;
                            total_snow_melt += snowmelt_m;
                            // Debug.WriteLine("GlobalMethods.t {0}, m {1} d {2} snowmelt {3} m", GlobalMethods.t,mo, day,snowmelt_m);
                            dailyflow(Pm[day], Dm[day], day, mo, snowmelt_m); // all snowmelt (+extra GlobalMethods.rain) becomes overland flow
                            P_OF += Pm[day];

                        }
                        else // no snow cover, rainfall intensity is used as threshold
                        {
                            if (Pm[day] > (Ks_min_mh * Dm[day]) && Dm[day] != 0) // develop. If rainfall is also spatially variable, this has to be adjusted
                            {
                                // Debug.WriteLine("wb4a");
                                //Debug.WriteLine("Overland flow initiated at date {0}/{1}/{2}", day, mo, GlobalMethods.t);
                                dailyflow(Pm[day], Dm[day], day, mo, 0);
                                P_OF += Pm[day];
                            }
                        }

                        if (snow_m < 0)
                        {
                            Debug.WriteLine("err_sno1");
                        }
                    }
                }

                //if (Pm.Sum() < P_OF) { MessageBox.Show("Pd > P"); }
                // Monthly water balance
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (GlobalMethods.dtm[row, col] != -9999)
                        {
                            double ET0m_act = ET0_m[mo] * Ra_rcm[row, col, mo] * veg_correction_factor[row, col];
                            // if (row == 0 & col == 0) { Debug.WriteLine(ET0m_act); }
                            ET0y[row, col] += ET0m_act;
                            double ETam = Pm.Sum() / Math.Pow(1 + Math.Pow(Pm.Sum() / (ET0m_act), 1.5), (1 / 1.5));

                            ETay[row, col] += ETam;

                            Iy[row, col] += (Pm.Sum() - P_OF) - ETam; // overland flow has been dealt with earlier (dailyflow), just as ponding
                            if (double.IsNaN(Iy[row, col]))
                            {
                                Debug.WriteLine("err wb1");
                            }
                        }
                    }
                }
                // Debugger.Break();

                daycount += month[mo];
            } // end months

            // yearly update infiltration

            bool Ineg = false;
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    if (Iy[row, col] < 0)
                    {
                        Ineg = true;
                    }
                }
            }
            if (Ineg)
            {
                ;
            }
            //if (GlobalMethods.t % 10 == 0) { Debugger.Break(); }
            // Debug.WriteLine("wb6");
        }

        void print_water_balance()
        {
            double P_wb = 0, ETa_wb = 0, I_wb = 0, OutF_wb = 0, snow_wb = 0, snow_left = 0, snow_start = 0, OFy_0 = 0, OFy_9 = 0;
            for (int rowb = 0; rowb < GlobalMethods.nr; rowb++)
            {
                for (int colb = 0; colb < GlobalMethods.nc; colb++)
                {
                    if (GlobalMethods.dtm[rowb, colb] != -9999)
                    {
                        P_wb += Py.Sum();
                        ETa_wb += ETay[rowb, colb];
                        I_wb += Iy[rowb, colb];
                        OutF_wb += total_outflow_y[rowb, colb];
                        snow_wb += snowfall_m;
                        snow_start += snow_start_m;
                        snow_left += snow_m;
                        OFy_0 += guiVariables.OFy_m[GlobalMethods.row, GlobalMethods.col, 0];
                        OFy_9 += guiVariables.OFy_m[GlobalMethods.row, GlobalMethods.col, 9];
                    }

                }
            }
            double balance = P_wb - ETa_wb - I_wb - OutF_wb;
            Debug.WriteLine("Snow start: {0}, snow end: {1}, P during snow: {2}", snow_start, snow_left, Pduringsnowmelt * GlobalMethods.nr * GlobalMethods.nc);
            Debug.WriteLine("Annual water balance");
            Debug.WriteLine("P: {0}, of which snowfall: {1}.  ETa: {2}. I: {3}. Outflow: {4}. Balance: {5}. OFy_0: {6}, OFy_9: {7}", P_wb, snow_wb, ETa_wb, I_wb, OutF_wb, balance, OFy_0, OFy_9);
            if (double.IsNaN(I_wb))
            {
                Debug.WriteLine("err_pwb1");
            }
            if (Math.Abs(balance + snow_start - snow_left) > 0.00000001)
            {
                Debug.WriteLine("err_pwb2");
            }
        }

        void print_spatial_water_balance()
        {
            for (int rowb = 0; rowb < GlobalMethods.nr; rowb++)
            {
                for (int colb = 0; colb < GlobalMethods.nc; colb++)
                {
                    Debug.WriteLine("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12}", rowb, colb, GlobalMethods.t, Py.Sum(), ETay[rowb, colb], guiVariables.OFy_m[rowb, colb, 0] - guiVariables.OFy_m[rowb, colb, 9], Iy[rowb, colb], (Py.Sum() - ETay[rowb, colb] - Iy[rowb, colb]), snow_m, snow_start_m, snowfall_m, total_outflow_y[rowb, colb], Pduringsnowmelt);
                    // snowfall_m is part of Py.Sum(), therefore, this should not be accounted for in the balance. Important terms are P, ETa, I and outflow. They close the balance. 
                }
            }
        }

        void print_P_ET0()
        {
            double ET_out = 0;
            double count = 0;
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    ET_out += ET0y[row, col];
                    count += 1;

                }
            }
            Debug.WriteLine("P {0} ET0 {1}", Py.Sum(), ET_out / count);
        }

        void update_potential_ET()
        {
            for (int hrow = 0; hrow < GlobalMethods.nr; hrow++)
            {
                for (int hcol = 0; hcol < GlobalMethods.nc; hcol++)
                {

                }
            }
        }

        void update_Ks() // both Ks matrix as Ks for topsoil
        {
            double[] tex_topsoil;
            double depth, fsilt, fclay, fOM, BD_t, slope_rad;
            int lay, topsoil;
            Ks_min_mh = 1000; Ks_max_mh = 0;
            List<double> BD_topsoil;


            // Debug.WriteLine("uks1");
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    BD_topsoil = new List<double>();
                    depth = 0;
                    tex_topsoil = new double[7];
                    lay = 0;
                    if (total_soil_mass(row, col) > 0)
                    {
                        while (depth <= 0.5 & lay < GlobalMethods.max_soil_layers)
                        {
                            // if (lay == GlobalMethods.max_soil_layers) { Debugger.Break(); }
                            if (total_layer_mass(row, col, lay) > 0)
                            {
                                depth += GlobalMethods.layerthickness_m[row, col, lay] / 2;

                                for (int text = 0; text < 5; text++)
                                {
                                    tex_topsoil[text] += GlobalMethods.texture_kg[row, col, lay, text];
                                }
                                tex_topsoil[5] += GlobalMethods.young_SOM_kg[row, col, lay];
                                tex_topsoil[5] += GlobalMethods.old_SOM_kg[row, col, lay];

                                BD_topsoil.Add(bulk_density_calc(GlobalMethods.texture_kg[row, col, lay, 0], GlobalMethods.texture_kg[row, col, lay, 1], GlobalMethods.texture_kg[row, col, lay, 2], GlobalMethods.texture_kg[row, col, lay, 3], GlobalMethods.texture_kg[row, col, lay, 4], GlobalMethods.old_SOM_kg[row, col, lay], GlobalMethods.young_SOM_kg[row, col, lay], depth));
                                depth += GlobalMethods.layerthickness_m[row, col, lay] / 2;

                            }
                            lay += 1;
                        }

                        fsilt = 100 * tex_topsoil[2] / (tex_topsoil[1] + tex_topsoil[2] + tex_topsoil[3] + tex_topsoil[4]); // only fine fraction
                        fclay = 100 * (tex_topsoil[3] + tex_topsoil[4]) / (tex_topsoil[1] + tex_topsoil[2] + tex_topsoil[3] + tex_topsoil[4]); // only fine fraction
                        fOM = 100 * tex_topsoil[5] / (tex_topsoil[1] + tex_topsoil[2] + tex_topsoil[3] + tex_topsoil[4] + tex_topsoil[5]); // only fine fraction
                        BD_t = BD_topsoil.Average() / 1000;
                        slope_rad = GlobalMethods.calc_slope_stdesc(GlobalMethods.row, GlobalMethods.col);
                        double slope_test = Math.Cos(slope_rad);

                        // Debug.WriteLine("uks3a");
                        Ks_topsoil_mh[row, col] = (Ks_wosten(fsilt, fclay, fOM, BD_t, 1) / 24) * Math.Cos(slope_rad);
                    }
                    else
                    {
                        Ks_topsoil_mh[row, col] = 0;
                        Debug.WriteLine("Empty soil at row {0}, col {1}, GlobalMethods.t {2}", row, col, GlobalMethods.t);
                    }

                    if (double.IsNaN(Ks_topsoil_mh[row, col]))
                    {
                        Debug.WriteLine("err_kst1");
                    }
                    // Debug.WriteLine("uks4");


                    //  Ks_topsoil_mh[row, col] *= veg_lag_factor;
                    // Ks_update
                    if (Ks_min_mh > Ks_topsoil_mh[row, col]) { Ks_min_mh = Ks_topsoil_mh[row, col]; }
                    if (Ks_max_mh < Ks_topsoil_mh[row, col]) { Ks_max_mh = Ks_topsoil_mh[row, col]; }
                    // Debug.WriteLine("uks5");
                }
            }
            ;
        }

        void create_daily_weather()
        {
            //Debug.WriteLine("dw.start");
            // 1. Select yearly timeseries and the corrected temperatures

            double P_ann = rain_value_m;
            int T_ann = temp_value_C;

            Random year_w = new Random(GlobalMethods.t); // GlobalMethods.t as random seed to get deterministic results
            int n_timeseries = year_w.Next(0, System.Convert.ToInt32(daily_n.Text));
            // n_timeseries = 5;

            //Debug.WriteLine("dw.1c");
            Array.Copy(guiVariables.Tmin_all, 365 * (n_timeseries), Tminy, 0, 365); // read new Tmin timeseries
            Array.Copy(guiVariables.Tmax_all, 365 * (n_timeseries), Tmaxy, 0, 365); // read new Tmin timeseries
            Array.Copy(guiVariables.Tavg_all, 365 * (n_timeseries), Tavgy, 0, 365); // read new Tmin timeseries

            Array.Copy(guiVariables.P_all, 365 * (n_timeseries), Py, 0, 365); // read new P timeseries
            Array.Copy(guiVariables.D_all, 365 * (n_timeseries), Dy, 0, 365); // read new D timeseries
            for (int pi = 0; pi < 365; pi++)
            {
                Py[pi] /= 1000; // convert to meters
                                // if (Py[pi] > 0.036) { Debugger.Break(); }
            }

            // 2. Rescale rainfall and temperature
            //Debug.WriteLine("dw.1a");

            if (check_scaling_daily_weather.Checked) // scaling with yearly values, for global change scenarios
            {
                double Py_sum = Py.Sum();
                for (int pi = 0; pi < Py.Count(); pi++) { Py[pi] = Py[pi] / Py_sum * P_ann; } // Scale with yearly P in meters
                double total_P = Py.Sum();

                int d_T = T_ann - Convert.ToInt32(Tavgy.Average());
                for (int pi = 0; pi < Tavgy.Count(); pi++)
                {
                    Tminy[pi] += d_T;
                    Tmaxy[pi] += d_T;
                    Tavgy[pi] += d_T;
                }
            }


            // 3. Calculate PET according to Hargreaves https://www.repository.utl.pt/bitstream/10400.5/4250/1/REP-J.L.Teixeira-InTech-Hargreaves_and_other_reduced_set_methods_for_calculating_evapotranspiration.pdf
            // Use monthly T values, gives a result very similar to daily values
            // multiplication with extraterrestrial radiation occurs in a later step, when ET0_m is actually necessary (water balance). Here we capture the monthly variation. At the later step, the spatiotemporal differences in solar radiation are captured
            int[] Tminm, Tmaxm, Tavgm;
            int daycount = 0;
            for (int pi = 0; pi < 12; pi++)
            {
                Tminm = new int[month[pi]];
                Tmaxm = new int[month[pi]];
                Tavgm = new int[month[pi]];

                Array.Copy(Tminy, daycount, Tminm, 0, month[pi]);
                Array.Copy(Tmaxy, daycount, Tmaxm, 0, month[pi]);
                Array.Copy(Tavgy, daycount, Tavgm, 0, month[pi]);

                daycount += month[pi];
                //get temperature info
                ET0_m[pi] = 0.0023 * (Tavgm.Average() + 17.78) * Math.Sqrt(Tmaxm.Average() - Tminm.Average()) * month[pi]; // multiply value by number of days in the month, to get the monthly total
                if (ET0_m[pi] < 0) { ET0_m[pi] = 0; }

            }


            //Debug.WriteLine("dw.end"); 
        }

        void dailyflow(double P_total, double D_total, int qday, int qmonth, double snowmelt)
        {
            // Debug.WriteLine("df1");
            /*
            Check:
            -is all water flow reset after first iteration?
            -do GlobalMethods.dtm and ponding values correspond well?
            -is the water balance closing?
            -do all cells refer to OFd[r,c,0] as flow component? rainfall and overland flow
            -link infiltration to differences in Ks (negative Ks) !!!
             */


            pond_d = new double[GlobalMethods.nr, GlobalMethods.nc];
            double[,] currentflow = new double[GlobalMethods.nr, GlobalMethods.nc];
            // Debug.WriteLine("df1");
            //every cell, inflow and outflow;
            double Qi, powered_slope_sum, OF_tot1 = 0, OF_tot2 = 0, OF_tot3 = 0; ;
            double[,,] OFd = new double[GlobalMethods.nr, GlobalMethods.nc, 10];

            //0: total flow
            //1-8: flow to neighbours
            //
            // 1 2 3 
            // 4   5
            // 6 7 8
            //
            //9: temporary flow, for the thresholds


            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    pond_d[row, col] = 0;
                    for (int it = 0; it <= 9; it++)
                    {
                        OFd[row, col, it] = 0;
                    } // reset all values
                }
            }


            int runner = 0;
            // Debug.WriteLine("df2");
            double totalwater = 0, totalwater2 = 0;
            // create overland flow in current flow map. This will be reset after flowing out, to consider a new flux of water after saddle overflow, without counting the first flux twice. 

            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    // Add rainfall excess to every cell. If negative, it can absorb incoming water from upstream.Overland flow will only be calculated when outflow is larger than zero. 
                    // Infiltration is dealt with at the end of the day. Water can still flow into the cell from higher up
                    //If there is snowmelt, all water (including rainfall), becomes overland flow
                    if (snowmelt > 0)
                    {
                        currentflow[row, col] = P_total + snowmelt;
                        total_water_flow += P_total + snowmelt;
                    }
                    else
                    {
                        currentflow[row, col] = P_total - Ks_topsoil_mh[row, col] * D_total; // infiltration excess becomes overland flow
                        total_water_flow += P_total - Ks_topsoil_mh[row, col] * D_total;
                    }
                    totalwater += currentflow[row, col];
                }
            }

            List<double> dh_list = new List<double>();
            List<string> dh_list_loc = new List<string>();
            // Debug.WriteLine("df3");
            // Debug.WriteLine("df2");
            // route to neighbours
            for (runner = GlobalMethods.number_of_data_cells - 1; runner >= 0; runner--)
            {
                //Debug.WriteLine("runner start of run: " + runner);
                if (GlobalMethods.index[runner] != -9999)
                {
                    int row = GlobalMethods.row_index[runner]; 
                    int col = GlobalMethods.col_index[runner];

                    //if (row == 186 & col == 499 & GlobalMethods.t == 5) { minimaps(186, 499); }


                    powered_slope_sum = 0;

                    // dh_list = new List<double>();

                    if (currentflow[row, col] > 0) // if there is currently water flowing out of the cell
                    {
                        for (int i = (-1); i <= 1; i++)
                        {
                            for (int j = (-1); j <= 1; j++)
                            {
                                dh = 0; dhtemp = -99999.99; GlobalMethods.d_x = GlobalMethods.dx;
                                if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                                {
                                    if (GlobalMethods.dtm[row + i, col + j] != -9999)
                                    {  //if the cell has no NODATA

                                        dh = GlobalMethods.dtm[row, col] - (GlobalMethods.dtm[row + i, col + j] + pond_d[row + i, col + j]);
                                        //Debug.WriteLine("dh = {0}", dh);
                                        if (dh > 0)
                                        {
                                            if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }   // for non-cardinal neighbours, we use the adapted length

                                            dh = dh / GlobalMethods.d_x;
                                            dh = Math.Pow(dh, conv_fac);
                                            dh_list.Add(dh);
                                            dh_list_loc.Add(Convert.ToString(col) + "." + Convert.ToString(j));
                                            powered_slope_sum = powered_slope_sum + dh;

                                            // no correction for possible sedimentation, like in normal overland flow
                                        }
                                    }
                                }
                            }
                        }
                        // Debug.WriteLine("df3");
                        if (powered_slope_sum > 0) // not in a depression
                        {
                            try
                            {

                                int direction = 0;
                                for (i = (-1); i <= 1; i++)
                                {
                                    for (j = (-1); j <= 1; j++)
                                    {
                                        if (!((i == 0) && (j == 0))) { direction++; }
                                        if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                                        {
                                            if (GlobalMethods.dtm[row + i, col + j] != -9999)
                                            {  //if the cell has no NODATA

                                                dh = GlobalMethods.dtm[row, col] - (GlobalMethods.dtm[row + i, col + j] + pond_d[row + i, col + j]);
                                                if (dh > 0)
                                                {
                                                    if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }   // for non-cardinal neighbours, we use the adapted length

                                                    dh = dh / GlobalMethods.d_x;
                                                    dh = Math.Pow(dh, conv_fac);
                                                    int t2 = GlobalMethods.t;
                                                    // if (col == 12) { Debugger.Break(); }
                                                    // flow to ij = dh / powered_slope_sum
                                                    OFd[row, col, direction] += dh / powered_slope_sum * currentflow[row, col]; // update outflow to each cell
                                                    currentflow[row + i, col + j] += dh / powered_slope_sum * currentflow[row, col];// update inflow of receiving cell. Here, a negative currentflow can become less negative or positive. 
                                                    OFd[row + i, col + j, 9] += dh / powered_slope_sum * currentflow[row, col]; // to track total inflow for the water balance

                                                }
                                            }
                                        }

                                    } // end j
                                } // end i
                                OFd[row, col, 0] += currentflow[row, col]; // total outflow
                                currentflow[row, col] = 0; // reset currentflow for possible later new flux after saddle position overflow
                            }
                            catch
                            {
                                MessageBox.Show("Error in water redistribution");
                            }


                        } // end powered_slope_sum > 0



                        else // no outflow, so at the edge of the catchment, or in a sink or depression
                        {
                            bool nodataneighbour = search_nodataneighbour(row, col);



                            if (nodataneighbour == true) // outflow
                            {
                                OFd[row, col, 0] += currentflow[row, col];
                                total_outflow_y[row, col] += currentflow[row, col];
                                currentflow[row, col] = 0;
                            }
                            else // depression
                            {
                                // Debug.WriteLine("df3a");
                                List<double> saddle_OF = new List<double>();
                                saddle_OF = ponding(row, col, currentflow[row, col], qday, qmonth);
                                OFd[row, col, 0] += currentflow[row, col];
                                currentflow[row, col] = 0;
                                // Debug.WriteLine("df3b");
                                // Debug.WriteLine("col =" + col);
                                if (saddle_OF.Count() > 0) // if there is saddle overflow
                                {
                                    int OF_row = Convert.ToInt32(saddle_OF[1]);
                                    int OF_col = Convert.ToInt32(saddle_OF[2]);
                                    currentflow[OF_row, OF_col] += saddle_OF[0];
                                    // Debug.WriteLine("runner old: " + runner + ", " + saddle_OF[0] + " " + saddle_OF[1] + " " + saddle_OF[2]);
                                    string rowcol = OF_row.ToString() + "." + OF_col.ToString();
                                    runner = Array.IndexOf(GlobalMethods.rowcol_index, rowcol) + 1; // GlobalMethods.index of selected row and col + 1, because in the next iteration of the for loop, 1 is subtracted from runner
                                                                                      // Debug.WriteLine("runner new: " + runner);
                                } // end of saddle overflow
                            }
                        }

                        // Debug.WriteLine("df4");
                    } // end inflow > 0
                    else
                    {


                        // two options:
                        // currentflow is still negative, indicating that only infiltration occurs. 
                        // currentflow has been reset earlier, so nothing happens
                        // if (snowmelt > 0) { Debugger.Break(); }
                        OFd[row, col, 0] += currentflow[row, col]; //Creates negative flow, because this is used for calculating infiltration. At that stage, this parameters will be set to zero
                        currentflow[row, col] = 0;
                    }


                } // end runner !=-9999
            } // end runner
              // Debug.WriteLine("df4");

            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    Iy[row, col] += pond_d[row, col];
                    if (double.IsNaN(Iy[row, col]))
                    {
                        Debug.WriteLine("err_df1");
                    }

                    // Calculate infiltration from rainfall and incoming flow
                    if (snowmelt == 0) // no snowmelt
                    {
                        if (OFd[row, col, 0] > 0) // if runoff occurred out of this cell, infiltration is maximum possible infiltration
                        {
                            Iy[row, col] += Ks_topsoil_mh[row, col] * D_total;
                            if (double.IsNaN(Iy[row, col]))
                            {
                                Debug.WriteLine("err_df2");
                            }

                        }
                        else
                        { // cell not saturated, infiltration is maximum infiltration minus the deficit at row, col, which is recorded in OFd[,,0]
                            Iy[row, col] += Ks_topsoil_mh[row, col] * D_total + OFd[row, col, 0]; //+, because OFd[row, col is negative
                                                                                                  //if (double.IsNaN(Iy[row, col])) { Debugger.Break(); }

                            OFd[row, col, 0] = 0;
                        }
                    }
                    else
                    {
                        // snow melt, no infiltration, because the soil is frozen. The rest leaves the catchment. this water only infiltrates when ponding. For the rest, it influences erosion
                    }



                    // Write daily flow to yearly flow
                    for (int it = 0; it <= 9; it++)
                    {
                        guiVariables.OFy_m[row, col, it] += OFd[row, col, it];
                    }

                    // check if total ponding equals total overland flow
                    pond_y[row, col] += pond_d[row, col];
                }
            }
            // Debug.WriteLine("df6");
        } // end dailyflow

        double Ks_wosten(double silt, double clay, double OM, double BD, int topsoil)
        {
            // Debug.WriteLine("KsW1");
            if (OM < 0.5) { OM = 0.5; } // half percent of OM for soils where it is absent, otherwise the PTFs will crash
            double KsW = Math.Exp(7.755 + 0.03252 * silt + 0.93 * topsoil - 0.967 * BD * BD - 0.000484 * clay * clay - 0.000322 * silt *
                silt + 0.001 / silt - 0.0748 / OM - 0.643 * Math.Log(silt) - 0.01398 * BD * clay - 0.1673 * BD * OM + 0.02986 * topsoil * clay - 0.03305 * topsoil * silt) / 100;
            // Debug.WriteLine("KsW2");
            if (Double.IsNaN(KsW)) { KsW = 0; }
            // Debug.WriteLine("KsW3");
            return (KsW); // m day-1

        }

        List<double> ponding(int row1, int col1, double inflow_m, int pday, int pmonth)
        {
            ponding_start = DateTime.Now;
            pond_d[row1, col1] += inflow_m;
            bool flatwater = false;
            int rowp, colp, rowp1, colp1, pondedcells, ri_of, ci_of;
            double minponding = inflow_m, dz_water;
            List<double> elev_p = new List<double>();
            List<double> output = new List<double>();
            elev_p.Add(GlobalMethods.dtm[row1, col1]);

            // 1. initiate list of ponded rows and cols
            List<int> pondingrows = new List<int>();
            List<int> pondingcols = new List<int>();
            pondingrows.Add(row1);
            pondingcols.Add(col1);
            List<double> dh_nb;
            // Debug.WriteLine("po1");

            // 2. loop over neighbours of ponded sink, to look for lowest neighbour to share water with
            // start with lowest neighbour of ponded cells, and work the way up
            while (flatwater == false)
            {
                dh_nb = new List<double>(); // store the dhs in a list, because we don'GlobalMethods.t know how big the pond gets and therefore how many neighbours there are

                for (int pondcell = 0; pondcell < pondingrows.Count(); pondcell++)
                {

                    rowp = pondingrows[pondcell];
                    colp = pondingcols[pondcell];
                    // if (minponding > pond_d[row, colp]) {minponding = pond_d[row, colp]; }
                    for (int i = -1; i <= 1; i++)
                    { // find lowest neighbour of all ponded cells.
                        for (int j = -1; j <= 1; j++)
                        {
                            if (((rowp + i) >= 0) && ((rowp + i) < GlobalMethods.nr) && ((colp + j) >= 0) && ((colp + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                            {
                                if (GlobalMethods.dtm[rowp + i, colp + j] != -9999)
                                {
                                    dh = (GlobalMethods.dtm[rowp, colp] + pond_d[rowp, colp]) - (GlobalMethods.dtm[rowp + i, colp + j] + pond_d[rowp + i, colp + j]);
                                    dh_nb.Add(dh);
                                } // end GlobalMethods.dtm != -9999
                            } // end if
                        } // end j
                    }  // end i
                }// end rowp colp

                //  Debug.WriteLine("po2");
                double maxdiff = dh_nb.Max();

                // 3. If there is a lower neighbour, share the water with him. Threshold a bit above 0, to avoid rounding errors
                if (maxdiff > 0.000000000001)
                {
                    // Debug.WriteLine("po3a");
                    for (int pondcell = 0; pondcell < pondingrows.Count(); pondcell++)
                    {
                        rowp = pondingrows[pondcell];
                        colp = pondingcols[pondcell];
                        double h_c = GlobalMethods.dtm[rowp, colp];
                        for (int i = -1; i <= 1; i++)
                        { // 4. relocate lowest neighbour of all ponded cells.
                            for (int j = -1; j <= 1; j++)
                            {
                                if (((rowp + i) >= 0) && ((rowp + i) < GlobalMethods.nr) && ((colp + j) >= 0) && ((colp + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                                {
                                    if (GlobalMethods.dtm[rowp + i, colp + j] != -9999)
                                    {
                                        // Debug.WriteLine("po3b");
                                        dh = (GlobalMethods.dtm[rowp, colp] + pond_d[rowp, colp]) - (GlobalMethods.dtm[rowp + i, colp + j] + pond_d[rowp + i, colp + j]);
                                        if (dh == dh_nb.Max()) // if the selected neighbour is the lowest neighbour
                                        {

                                            pondedcells = pondingrows.Count();
                                            dz_water = pondedcells * dh / (pondedcells + 1); // amount of water added to the lowest (waterless) neighbour

                                            // if change in ponded cells exceeds available water in shallowest ponding cell (e.g. saddle positions), water drop is limited to this depth, preventing negative ponding
                                            double h_n = GlobalMethods.dtm[rowp + i, colp + j] + pond_d[rowp + i, colp + j];
                                            if (h_n < h_c) // if we cross a saddle position, stop water flow for now

                                            {
                                                double overflow = pond_d[rowp, colp] * pondedcells; // excess water is ponding level at the saddle position, times the amount of ponding cells

                                                for (int i_of = 0; i_of < pondedcells; i_of++)
                                                {
                                                    ri_of = pondingrows[i_of];
                                                    ci_of = pondingcols[i_of];
                                                    pond_d[ri_of, ci_of] -= overflow / pondedcells;
                                                }
                                                output.Add(overflow); // create list for output, consisting of saddle r and c, and amount of overflow
                                                output.Add(rowp);
                                                output.Add(colp);
                                                //OFd[rowp, colp, 0] += overflow; // add overflow to the saddle cell, not the lower neighbour. from the saddle cell it will flow onward
                                                dz_water = 0; flatwater = true;
                                            }

                                            else
                                            {
                                                // Debug.WriteLine("po3c");
                                                for (int pondcell1 = 0; pondcell1 < pondingrows.Count(); pondcell1++)
                                                { // subtract redistributed water from original ponds
                                                    rowp1 = pondingrows[pondcell1];
                                                    colp1 = pondingcols[pondcell1];
                                                    pond_d[rowp1, colp1] -= dz_water / pondedcells;
                                                }
                                                pond_d[rowp + i, colp + j] += dz_water; // add water to the new ponding cell
                                                                                        //  Debug.WriteLine("po3d");
                                                pondingrows.Add(rowp + i);
                                                pondingcols.Add(colp + j);
                                                if (minponding > pond_d[rowp + i, colp + j]) { minponding = pond_d[rowp + i, colp + j]; }
                                                elev_p.Add(GlobalMethods.dtm[rowp + i, colp + j]);
                                            }


                                            //for(int iii = 0;iii<pondingrows.Count();iii++)
                                            //{
                                            //    Debug.Write(GlobalMethods.dtm[pondingrows[iii], pondingcols[iii]] + pond_d[pondingrows[iii], pondingcols[iii]]+" ");
                                            //}
                                            //Debug.Write("\n");
                                        }
                                    } // end GlobalMethods.dtm != -9999
                                } // end if
                            } // end j
                        }  // end i

                    } // end rowp colp
                } // end if dh_nb > 
                else
                { // no lower neighbour, end of ponding
                    flatwater = true;
                }
                // Debug.WriteLine("po3");
            } // end of water redistribution


            bool negativeponding = false;
            for (int r = 0; r < GlobalMethods.nr; r++)
            {
                for (int c = 0; c < GlobalMethods.nc; c++)
                {
                    //pond_y[r, c] += pond_d[r, c];
                    if (pond_d[r, c] < -0.00000001) { Debug.WriteLine("negative ponding on row {0} and col {1}. Amount = {2} ", row1, col1, pond_d[r, c]); }
                }
            }
            // if(negativeponding == true) { GlobalMethods.out_double(GlobalMethods.Workdir + "\\debug\\" + run_number + "_" + GlobalMethods.t + "_ "+ pmonth+"_ " + pday + "_out_ponding.asc", pond_d); }
            ponding_t += DateTime.Now - ponding_start;

            return (output);

        }

        bool stagnation(double I_d)
        {
            // DEVELOP MM account for different moisture conditions throughout the year
            bool stag;
            bool stag_total = false;
            int lay;
            double depth;
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    if (GlobalMethods.dtm[row, col] != -9999)
                    {
                        stag = false;
                        lay = 0;
                        depth = 0;
                        if (I_d > Ks_md[row, col, lay]) { stag = true; stag_total = true; }
                        while (stag == false)
                        {
                            depth += GlobalMethods.layerthickness_m[row, col, lay];
                            if (I_d > Ks_md[row, col, lay]) { stag = true; stag_total = true; }
                        }
                        if (stag == true) { stagdepth[row, col] = depth; }
                    }
                }
            }
            return (stag_total);
        }

        void lateral_flow()
        {

            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    for (i = (-1); i <= 1; i++)
                    {
                        for (j = (-1); j <= 1; j++)
                        {
                            if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                            {
                                if (GlobalMethods.dtm[row + i, col + j] != -9999)
                                {  //if the cell has no NODATA

                                }
                            }
                        }
                    }
                }
            }
        }

        #endregion

        #region Soil development code

        void update_all_soil_thicknesses(int row_update, int col_update)
        {
            for (int layer_update = 0; layer_update < GlobalMethods.max_soil_layers; layer_update++)
            {
                GlobalMethods.layerthickness_m[row_update, col_update, layer_update] = thickness_calc(row_update, col_update, layer_update);
                GlobalMethods.layerthickness_m[row_update, col_update, layer_update] = thickness_calc(row_update, col_update, layer_update);
            }

        }

        void remove_empty_layers(int row2, int col2)
        {
            // mainly after tree fall, there can be empty soil layers at the surface. This module shifts the layers up.
            // displaysoil(row2, col2);
            // Debug.WriteLine("rel1");
            // DEVELOP after shifting cells up the script runs through the lower empty cells, moving them up also. with some booleans, this should be prevented. there is no error now, only longer simulation time.  
            try
            {
                if (diagnostic_mode == 1) { Debug.WriteLine("entered removing empty layers"); }
                int empty_layers = 0;
                bool shift_layers = false;

                int n_shifts = 0;
                double mass_before = total_soil_mass(row2, col2);
                // Debug.WriteLine("rel2");
                for (int lay2 = 0; lay2 < GlobalMethods.max_soil_layers; lay2++)
                {
                    bool full_layer_shift = false;
                    if (total_layer_mass(row2, col2, lay2) < 0.000000000001) // empty layer
                                                                             // mind for empty layers at the bottom
                    {
                        shift_layers = true;
                        // if(n_shifts == 0) { displaysoil(row2, col2); }
                        // n_shifts += 1;

                        empty_layers++;
                        for (int layert = lay2; layert < GlobalMethods.max_soil_layers - 1; layert++) // for all underlying layers, shift one up (since there is space anyway)
                        {

                            if (total_layer_mass(row2, col2, layert + 1) > 0) { full_layer_shift = true; }
                            // Debug.WriteLine(layert);
                            for (i = 0; i < 5; i++)
                            {
                                GlobalMethods.texture_kg[row2, col2, layert, i] = GlobalMethods.texture_kg[row2, col2, layert + 1, i];
                            }
                            GlobalMethods.old_SOM_kg[row2, col2, layert] = GlobalMethods.old_SOM_kg[row2, col2, layert + 1];
                            GlobalMethods.young_SOM_kg[row2, col2, layert] = GlobalMethods.young_SOM_kg[row2, col2, layert + 1];
                        }
                        for (i = 0; i < 5; i++)
                        {
                            GlobalMethods.texture_kg[row2, col2, GlobalMethods.max_soil_layers - 1, i] = 0;
                        }
                        GlobalMethods.old_SOM_kg[row2, col2, GlobalMethods.max_soil_layers - 1] = 0;
                        GlobalMethods.young_SOM_kg[row2, col2, GlobalMethods.max_soil_layers - 1] = 0;
                        GlobalMethods.layerthickness_m[row2, col2, GlobalMethods.max_soil_layers - 1] = 0;

                        if (full_layer_shift == true) { lay2--; }
                        //Debug.WriteLine("-");
                        //displaysoil(row2, col2);
                    }
                }
                // Debug.WriteLine("rel3");
                if (shift_layers == true)
                {
                    double mass_after = total_soil_mass(row2, col2);
                    if (Math.Round(mass_before - mass_after) > 0.0000001)
                    {
                        Debug.WriteLine("Loss of soil data during removal of empty layers");
                        displaysoil(row2, col2);
                    }

                }

                if (diagnostic_mode == 1) { if (n_shifts > 0) { displaysoil(row2, col2); Debug.WriteLine("n layers shifted: {0} in row {1}, col {2}", n_shifts, row2, col2); } }
                update_all_soil_thicknesses(row2, col2);
            } // end try

            catch
            {
                Debug.WriteLine("Error in removing empty layers");
            }
            if (diagnostic_mode == 1) { Debug.WriteLine(" at end of remove_empty_layers: "); displaysoil(row2, col2); }
        }

        /*void soil_update_split_and_combine_layers()
        {

            //where at the end of soil development, splitting and combining of soil layers is performed.

            //per layer: if too thin: combine with one of the two neighbours (the closest one in properties). 
            // too thick: split
            // if total too much - combine the most similar two layers although the product conforms to most restrictive rule about thicknesses
            //maat voor verschil is som (absolute differences in de vijf text classes and two organic matter classes)

            // 0-50 cm    min 2.5   insteek 5    maximum 10 cm       n=10    bovenste laag geen minimum (sediment HOEFT niet meteen weggemiddeld te worden - pas als nodig)
            // 50-250 cm  min 10    insteek 25    maximum 50 cm      n=8
            // daarna     min 50    insteek 100  geen max            n=4

            //for combining : only when matrix size exceeded, then those two layers that are most equal are combined
            //do not combine two layers that together are too thick for their position in the profile. If needed, make lowest layer thicker. Therefore, approach from above

            //MISSING: HOW DO WE GO FROM MAX TO LESS LAYERS? ART
          //  Debug.WriteLine("suscl1 ");

            //displaysoil(0,0);

            total_average_soilthickness_m = 0;
            number_soil_thicker_than = 0;
            number_soil_coarser_than = 0;
            local_soil_depth_m = 0;
            local_soil_mass_kg = 0;

            int layer;
            int numberoflayers = 0;
            double depth_m;  // keep better track of this, currently not OK yet
            try
            {
                //Debug.WriteLine("suscl2");
                //displaysoil(0, 0);
                for (row = 0; row < GlobalMethods.nr; row++)
                {
                    for (col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (GlobalMethods.dtm[row, col] != -9999)
                        {
                            GlobalMethods.soildepth_m[row, col] = total_soil_thickness(row, col);
                            GlobalMethods.soildepth_m[row, col] = total_soil_thickness(row, col);

                           // if (GlobalMethods.soildepth_m[row, col] < 50) { Debug.WriteLine("susac d = {0}, at GlobalMethods.t {1}, r {2}, c {3}", GlobalMethods.soildepth_m[row, col], GlobalMethods.t, row, col); displaysoil(row, col); Debugger.Break(); }
                            depth_m = 0;
                            for (layer = 0; layer < GlobalMethods.max_soil_layers; layer++)
                            {
                                ////update the layers' thickness now that textures and organic matter amounts have changed (if there is anything in the layer at all).
                                //if (!(GlobalMethods.texture_kg[row, col, layer, 0] == 0 && GlobalMethods.texture_kg[row, col, layer, 1] == 0 && GlobalMethods.texture_kg[row, col, layer, 2] == 0 && GlobalMethods.texture_kg[row, col, layer, 3] == 0 && GlobalMethods.texture_kg[row, col, layer, 4] == 0 && GlobalMethods.young_SOM_kg[row, col, layer] == 0 && GlobalMethods.old_SOM_kg[row, col, layer] == 0))
                                //{
                                GlobalMethods.layerthickness_m[row, col, layer] = thickness_calc(row, col, layer);
                                //}
                                if (timeseries.timeseries_soil_mass_checkbox.Checked && System.Convert.ToInt32(timeseries.timeseries_soil_cell_row.Text) == row && System.Convert.ToInt32(timeseries.timeseries_soil_cell_col.Text) == col)
                                {
                                    local_soil_mass_kg += GlobalMethods.texture_kg[row, col, layer, 0] + GlobalMethods.texture_kg[row, col, layer, 1] + GlobalMethods.texture_kg[row, col, layer, 2] + GlobalMethods.texture_kg[row, col, layer, 3] + GlobalMethods.texture_kg[row, col, layer, 4] + GlobalMethods.young_SOM_kg[row, col, layer] + GlobalMethods.old_SOM_kg[row, col, layer];
                                }
                                if (timeseries.timeseries_soil_mass_checkbox.Checked && layer == 0 && GlobalMethods.texture_kg[row, col, layer, 0] / (GlobalMethods.texture_kg[row, col, layer, 0] + GlobalMethods.texture_kg[row, col, layer, 1] + GlobalMethods.texture_kg[row, col, layer, 2] + GlobalMethods.texture_kg[row, col, layer, 3] + GlobalMethods.texture_kg[row, col, layer, 4] + GlobalMethods.young_SOM_kg[row, col, layer] + GlobalMethods.old_SOM_kg[row, col, layer]) > System.Convert.ToDouble(timeseries.timeseries_soil_coarser_fraction_textbox.Text))
                                {
                                    number_soil_coarser_than++;
                                }
                            }
                        }

                    }
                }


                //displaysoil(0, 0);
                for (row = 0; row < GlobalMethods.nr; row++)
                {
                    for (col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (GlobalMethods.dtm[row, col] != -9999)
                        {
                            // Debug.WriteLine("suscl0" + row + ", " + col + ", " + GlobalMethods.t + " " + total_soil_mass(row, col));
                            //Debug.WriteLine("soil before splitting");
                            // Debug.WriteLine("uscl1");
                            // if (row == 0 & col == 0) { displaysoil(row, col); }
                            depth_m = 0; numberoflayers = 0;


                            double mass1 = total_soil_mass(row, col);
                            for (layer = 0; layer < GlobalMethods.max_soil_layers; layer++)
                            {
                                if (GlobalMethods.layerthickness_m[row, col, layer] > 0)
                                {
                                    depth_m += GlobalMethods.layerthickness_m[row, col, layer];
                                    //Debug.WriteLine("depth is now " + depth + " for lyr " +  layer);
                                    numberoflayers++;

                                    // 0-50 cm    min 2.5   insteek 5    maximum 10 cm       n=10    bovenste laag geen minimum (sediment HOEFT niet meteen weggemiddeld te worden - pas als nodig)
                                    // 50-250 cm  min 10    insteek 25    maximum 50 cm      n=8
                                    // daarna     min 50    insteek 100  geen max            n=4

                                    if (depth_m <= 0.5)
                                    {

                                        if (GlobalMethods.layerthickness_m[row, col, layer] < 0.025 && layer != 0)
                                        { //combine layers: select the one most like this one
                                            if (layer_difference(row, col, layer, layer - 1) > layer_difference(row, col, layer, layer + 1))
                                            {
                                                combine_layers(row, col, layer, layer + 1);
                                                update_all_soil_thicknesses(row, col);
                                                numberoflayers--;

                                            }
                                            else
                                            {
                                                combine_layers(row, col, layer - 1, layer);
                                                layer--;  //because we combined with the previous one, the current one has been replaced with one that has not yet been considered
                                                update_all_soil_thicknesses(row, col);
                                                numberoflayers--;

                                            }
                                            // we will now check whether layers have become too thick and if needed cut them in half
                                        }

                                        if (GlobalMethods.layerthickness_m[row, col, layer] > 0.1)
                                        { //split 
                                          // Debug.WriteLine("splitting after combining 1");
                                            split_layer(row, col, layer, depth_m);
                                            update_all_soil_thicknesses(row, col);
                                            numberoflayers++;

                                        }

                                        // 0-50 cm    min 2.5   insteek 5    maximum 10 cm       n=10    bovenste laag geen minimum (sediment HOEFT niet meteen weggemiddeld te worden - pas als nodig)
                                        // 50-250 cm  min 10    insteek 25    maximum 50 cm      n=8
                                        // daarna     min 50    insteek 100  geen max            n=4

                                    }

                                    // Debug.WriteLine("uscl2");

                                    if (depth_m > 0.5 && depth_m <= 3.0)
                                    {
                                        if (GlobalMethods.layerthickness_m[row, col, layer] < 0.1 && layer != 0)
                                        { //combine 
                                            if (layer_difference(row, col, layer, layer - 1) > layer_difference(row, col, layer, layer + 1))
                                            {
                                                combine_layers(row, col, layer, layer + 1);
                                                update_all_soil_thicknesses(row, col);
                                                numberoflayers--;

                                            }
                                            else
                                            {
                                                combine_layers(row, col, layer - 1, layer);
                                                layer--;  //because we combined with the previous one, the current one has been replaced with one that has not yet been considered
                                                update_all_soil_thicknesses(row, col);
                                            }
                                        }
                                        if (GlobalMethods.layerthickness_m[row, col, layer] > 0.5)
                                        { //split 
                                          // Debug.WriteLine("splitting after combining 2. layer = {0}, layerthickness = {1}", layer, GlobalMethods.layerthickness_m[row, col, layer]);
                                            split_layer(row, col, layer, depth_m);
                                            update_all_soil_thicknesses(row, col);
                                        }
                                    }

                                    // Debug.WriteLine("uscl3");

                                    if (depth_m > 3.0)
                                    {
                                        // Debug.WriteLine("uscl3");
                                        if (GlobalMethods.layerthickness_m[row, col, layer] < 0.5 && layer != 0 && layer<GlobalMethods.max_soil_layers)
                                        { //combine 
                                            Debug.WriteLine("uscl3");
                                            if (layer_difference(row, col, layer, layer - 1) > layer_difference(row, col, layer, layer + 1))
                                            {
                                                // Debug.WriteLine("uscl3.1");
                                                combine_layers(row, col, layer, layer + 1);
                                                update_all_soil_thicknesses(row, col);
                                                numberoflayers--;
                                                // Debug.WriteLine("uscl3.2");


                                            }
                                            else
                                            {
                                                Debug.WriteLine("uscl3.3");

                                                combine_layers(row, col, layer - 1, layer);
                                                layer--;  //because we combined with the previous one, the current one has been replaced with one that has not yet been considered
                                                numberoflayers--;
                                                update_all_soil_thicknesses(row, col);
                                                // Debug.WriteLine("uscl3.4");

                                            }

                                        }
                                        if (GlobalMethods.layerthickness_m[row, col, layer] > 0.5)
                                        { //split 
                                          // no splitting, no maximum thickness
                                        }
                                    }
                                    // Debug.WriteLine("uscl4");

                                    //Debug.WriteLine("depth is now " + depth + " and number of layers is  " + numberoflayers);
                                }
                            }

                            //Debug.WriteLine("suscl1" + row + ", " + col + ", " + GlobalMethods.t);
                            // Debug.WriteLine("uscl5");

                            for (int layupdate = 0; layupdate < GlobalMethods.max_soil_layers; layupdate++)
                            {
                                GlobalMethods.layerthickness_m[row, col, layupdate] = thickness_calc(row, col, layupdate);
                            }

                            double mass2 = total_soil_mass(row, col);
                            // Debug.WriteLine("uscl6");

                            //Debug.WriteLine("Soil after splitting");
                            //if (row == 0 & col == 0) { displaysoil(row, col); }
                            //Debug.WriteLine("suscl4");
                            //displaysoil(0, 0);
                            //if (numberoflayers > GlobalMethods.max_soil_layers)
                            if (Math.Abs(mass1 - mass2) > 0.0000001) // meij. changed the check to amount of material in the soil. with a lot of combining and splitting, the count of numberoflayers goes wrong. Threshold not 0, due to rounding errors
                            {
                               // this should never happen, because the data of the lowest layers have then been lost.
                               Debug.WriteLine(" Warning - loss of soil data "); 
                               Debug.WriteLine("mass difference = {0}", (mass1 - mass2)); 
                               displaysoil(row, col);
                               Debugger.Break();
                            }
                            // Debug.WriteLine("uscl7");

                            if (timeseries.timeseries_number_soil_thicker_checkbox.Checked && System.Convert.ToDouble(timeseries.timeseries_soil_thicker_textbox.Text) < depth_m) { number_soil_thicker_than++; }
                            if (timeseries.total_average_soilthickness_checkbox.Checked) { total_average_soilthickness_m += depth_m; }
                            if (timeseries.timeseries_soil_depth_checkbox.Checked && System.Convert.ToInt32(timeseries.timeseries_soil_cell_row.Text) == row && System.Convert.ToInt32(timeseries.timeseries_soil_cell_col.Text) == col)
                            {
                                local_soil_depth_m = depth_m;
                            }

                            //Debug.WriteLine("suscl0" + row + ", " + col + ", " + GlobalMethods.t + " " + total_soil_mass(row, col));
                            //double old_thickness = GlobalMethods.soildepth_m[row, col];
                            //double new_thickness = total_soil_thickness(row, col);
                            //double dtm_difference = GlobalMethods.soildepth_m[row, col] - total_soil_thickness(row, col);
                            //if (Math.Abs(old_thickness - new_thickness) > 0.0005)
                            //{
                            //    displaysoil(row, col);
                            //    Debugger.Break();
                            //}
                            //double old_thickness = GlobalMethods.soildepth_m[row, col];
                            //double new_thickness = total_soil_thickness(row, col);
                            //if(Math.Abs(old_thickness-new_thickness)>0.5)
                            //{
                            //    Debugger.Break();
                            //}
                            GlobalMethods.dtm[row, col] += GlobalMethods.soildepth_m[row, col] - total_soil_thickness(row, col); // update GlobalMethods.dtm by differences in soil depth
                            GlobalMethods.soildepth_m[row, col] = total_soil_thickness(row, col);

                        }

                    }
                }
               if (timeseries.total_average_soilthickness_checkbox.Checked) { total_average_soilthickness_m /= GlobalMethods.number_of_data_cells; }

            }
            catch
            {
                Debug.WriteLine("Error in updating, splitting and combining layers");
            }

        }
         * */ //aangepast voor tree fall, maar met fout erin

        void soil_update_split_and_combine_layers_standard()
        {
            //where at the end of soil development, splitting and combining of soil layers is performed.

            //per layer: if too thin: combine with one of the two neighbours (the closest one in properties). 
            // too thick: split
            // if total too much - combine the most similar two layers although the product conforms to most restrictive rule about thicknesses
            //maat voor verschil is som (absolute differences in de vijf text classes and two organic matter classes)

            // 0-50 cm    min 2.5   insteek 5    maximum 10 cm       n=10    bovenste laag geen minimum (sediment HOEFT niet meteen weggemiddeld te worden - pas als nodig)
            // 50-200 cm  min 10    insteek 15    maximum 50 cm      n=8
            // daarna     min 40    insteek 50  geen max            n=4

            //for combining : only when matrix size exceeded, then those two layers that are most equal are combined
            //do not combine two layers that together are too thick for their position in the profile. If needed, make lowest layer thicker. Therefore, approach from above

            //MISSING: HOW DO WE GO FROM MAX TO LESS LAYERS? ART
            //  Debug.WriteLine("suscl1 ");

            //displaysoil(0,0);
            double mass_before = total_catchment_mass();
            total_average_soilthickness_m = 0;
            number_soil_thicker_than = 0;
            number_soil_coarser_than = 0;
            local_soil_depth_m = 0;
            local_soil_mass_kg = 0;

            int layer;
            int numberoflayers = 0;
            double depth_m;  // keep better track of this, currently not OK yet
            try
            {
                //Debug.WriteLine("suscl2");
                //displaysoil(0, 0);
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (GlobalMethods.dtm[row, col] != -9999)
                        {
                            update_all_soil_thicknesses(row, col);
                            update_all_soil_thicknesses(row, col);

                            depth_m = 0;
                            for (layer = 0; layer < GlobalMethods.max_soil_layers; layer++)
                            {
                                ////update the layers' thickness now that textures and organic matter amounts have changed (if there is anything in the layer at all).
                                //if (!(GlobalMethods.texture_kg[row, col, layer, 0] == 0 && GlobalMethods.texture_kg[row, col, layer, 1] == 0 && GlobalMethods.texture_kg[row, col, layer, 2] == 0 && GlobalMethods.texture_kg[row, col, layer, 3] == 0 && GlobalMethods.texture_kg[row, col, layer, 4] == 0 && GlobalMethods.young_SOM_kg[row, col, layer] == 0 && GlobalMethods.old_SOM_kg[row, col, layer] == 0))
                                //{
                                //GlobalMethods.layerthickness_m[row, col, layer] = thickness_calc(row, col, layer);
                                //GlobalMethods.layerthickness_m[row, col, layer] = thickness_calc(row, col, layer);
                                //if (GlobalMethods.layerthickness_m[row,col,layer] < 0) { Debugger.Break(); } //MMS
                                //find_negative_texture_rcl(row, col, layer); //MMS
                                //}
                                if (guiVariables.Timeseries.Timeseries_soil_mass_checkbox && System.Convert.ToInt32(guiVariables.Timeseries.Timeseries_soil_cell_row) == row && System.Convert.ToInt32(guiVariables.Timeseries.Timeseries_soil_cell_col) == col)
                                {
                                    local_soil_mass_kg += GlobalMethods.texture_kg[row, col, layer, 0] + GlobalMethods.texture_kg[row, col, layer, 1] + GlobalMethods.texture_kg[row, col, layer, 2] + GlobalMethods.texture_kg[row, col, layer, 3] + GlobalMethods.texture_kg[row, col, layer, 4] + GlobalMethods.young_SOM_kg[row, col, layer] + GlobalMethods.old_SOM_kg[row, col, layer];
                                    if (local_soil_mass_kg < 0)
                                    {
                                        Debug.WriteLine("err_uscl1");
                                    } //MMS
                                }
                                if (guiVariables.Timeseries.Timeseries_soil_mass_checkbox && layer == 0 && GlobalMethods.texture_kg[row, col, layer, 0] / (GlobalMethods.texture_kg[row, col, layer, 0] + GlobalMethods.texture_kg[row, col, layer, 1] + GlobalMethods.texture_kg[row, col, layer, 2] + GlobalMethods.texture_kg[row, col, layer, 3] + GlobalMethods.texture_kg[row, col, layer, 4] + GlobalMethods.young_SOM_kg[row, col, layer] + GlobalMethods.old_SOM_kg[row, col, layer]) > System.Convert.ToDouble(guiVariables.Timeseries.Timeseries_soil_coarser_fraction_textbox))
                                {
                                    number_soil_coarser_than++;
                                }
                            }

                        }
                    }
                }


                //displaysoil(0, 0);
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (GlobalMethods.dtm[row, col] != -9999)
                        {
                            double old_soil_mass = total_soil_mass(row, col), new_soil_mass;
                            // Debug.WriteLine("suscl0" + row + ", " + col + ", " + GlobalMethods.t + " " + total_soil_mass(row, col));
                            //Debug.WriteLine("soil before splitting");
                            // if (row == 0 & col == 0) { displaysoil(row, col); }
                            depth_m = 0; numberoflayers = 0;
                            bool boolsplit = false;
                            bool boolcombine = false;
                            for (layer = 0; layer < GlobalMethods.max_soil_layers - 1; layer++)
                            {
                                if (GlobalMethods.layerthickness_m[row, col, layer] > 0)
                                {

                                    //Debug.WriteLine("depth is now " + depth + " for lyr " +  layer);
                                    numberoflayers++;

                                    // 0-50 cm    min 2.5   insteek 5    maximum 10 cm       n=10    bovenste laag geen minimum (sediment HOEFT niet meteen weggemiddeld te worden - pas als nodig)
                                    // 50-250 cm  min 10    insteek 25    maximum 50 cm      n=8
                                    // daarna     min 50    insteek 100  geen max            n=4

                                    if (depth_m <= 0.5)
                                    {
                                        if (layer == 0 & GlobalMethods.layerthickness_m[row, col, layer] < 0.001) // smaller than one mm -> merge with layer below
                                        {
                                            combine_layers(row, col, layer, layer + 1);
                                            update_all_soil_thicknesses(row, col);
                                            boolcombine = true;
                                            if (Math.Round(old_soil_mass, 6) != Math.Round(total_soil_mass(row, col), 6)) { Debug.WriteLine("err_uscl2"); }
                                        }

                                        if (GlobalMethods.layerthickness_m[row, col, layer] < 0.025 && layer != 0)
                                        { //combine layers: select the one most like this one
                                            if (layer_difference(row, col, layer, layer - 1) > layer_difference(row, col, layer, layer + 1))
                                            {
                                                combine_layers(row, col, layer, layer + 1);
                                                update_all_soil_thicknesses(row, col);
                                                boolcombine = true;
                                                if (Math.Round(old_soil_mass, 6) != Math.Round(total_soil_mass(row, col), 6))
                                                {
                                                    Debug.WriteLine("err_uscl3");
                                                }

                                            }
                                            else
                                            {
                                                combine_layers(row, col, layer - 1, layer);
                                                layer--;  //because we combined with the previous one, the current one has been replaced with one that has not yet been considered
                                                update_all_soil_thicknesses(row, col);
                                                boolcombine = true;
                                                // if (Math.Round(old_soil_mass, 6) != Math.Round(total_soil_mass(row, col), 6)) { Debugger.Break(); }

                                            }
                                            // we will now check whether layers have become too thick and if needed cut them in half
                                        }

                                        while (GlobalMethods.layerthickness_m[row, col, layer] > 0.1)
                                        { //split 
                                          // Debug.WriteLine("splitting after combining 1");
                                            split_layer(row, col, layer, depth_m);
                                            update_all_soil_thicknesses(row, col);
                                            boolsplit = true;
                                            // if (Math.Abs(old_soil_mass-total_soil_mass(row, col))>0.00000001) { Debugger.Break(); }

                                        }

                                        // 0-50 cm    min 2.5   insteek 5    maximum 10 cm       n=10    bovenste laag geen minimum (sediment HOEFT niet meteen weggemiddeld te worden - pas als nodig)
                                        // 50-250 cm  min 10    insteek 25    maximum 50 cm      n=8
                                        // daarna     min 50    insteek 100  geen max            n=4

                                    }

                                    depth_m += GlobalMethods.layerthickness_m[row, col, layer]; // MM moved this down one if-function, to be able to split big clumps of earth by tree fall. If I put it at the end, it will give problems with splitting the one-before-last layer




                                    if (depth_m > 0.5 && depth_m <= 2)
                                    {
                                        if (GlobalMethods.layerthickness_m[row, col, layer] < 0.1 && layer != 0)
                                        { //combine 
                                            if (layer_difference(row, col, layer, layer - 1) > layer_difference(row, col, layer, layer + 1))
                                            {
                                                combine_layers(row, col, layer, layer + 1);
                                                update_all_soil_thicknesses(row, col);
                                                boolcombine = true;
                                                if (Math.Abs(old_soil_mass - total_soil_mass(row, col)) > 0.00000001)
                                                {
                                                    Debug.WriteLine("err_uscl4");
                                                }
                                            }
                                            else
                                            {
                                                combine_layers(row, col, layer - 1, layer);
                                                layer--;  //because we combined with the previous one, the current one has been replaced with one that has not yet been considered
                                                update_all_soil_thicknesses(row, col);
                                                boolcombine = true;
                                                // if (Math.Abs(old_soil_mass - total_soil_mass(row, col)) > 0.00000001) { Debugger.Break(); }

                                            }
                                        }
                                        while (GlobalMethods.layerthickness_m[row, col, layer] > 0.5)
                                        { //split 
                                          // Debug.WriteLine("splitting after combining 2. layer = {0}, layerthickness = {1}", layer, GlobalMethods.layerthickness_m[row, col, layer]);
                                            split_layer(row, col, layer, depth_m);
                                            update_all_soil_thicknesses(row, col);
                                            boolsplit = true;
                                            new_soil_mass = total_soil_mass(row, col);
                                            if (Math.Abs(old_soil_mass - new_soil_mass) > 0.00000001)
                                            {
                                                Debug.WriteLine("err_uscl5");
                                            }

                                        }
                                    }


                                    if (depth_m > 2)
                                    {
                                        if (GlobalMethods.layerthickness_m[row, col, layer] < 0.4 && layer != 0)
                                        { //combine 
                                          //displaysoil(row, col);
                                            if (layer_difference(row, col, layer, layer - 1) > layer_difference(row, col, layer, layer + 1))
                                            {
                                                combine_layers(row, col, layer, layer + 1);
                                                update_all_soil_thicknesses(row, col);
                                                boolcombine = true;
                                                if (Math.Abs(old_soil_mass - total_soil_mass(row, col)) > 0.00000001)
                                                {
                                                    Debug.WriteLine("err_uscl6");
                                                }

                                            }
                                            else
                                            {
                                                combine_layers(row, col, layer - 1, layer);
                                                layer--;  //because we combined with the previous one, the current one has been replaced with one that has not yet been considered
                                                numberoflayers--;
                                                update_all_soil_thicknesses(row, col);
                                                boolcombine = true;
                                                // if (Math.Abs(old_soil_mass - total_soil_mass(row, col)) > 0.00000001) { Debugger.Break(); }

                                            }

                                        }
                                        if (GlobalMethods.layerthickness_m[row, col, layer] > 0.5)
                                        { //split 
                                          // no splitting, no maximum thickness
                                        }
                                    }


                                    //Debug.WriteLine("depth is now " + depth + " and number of layers is  " + numberoflayers);
                                }
                            }
                            //Debug.WriteLine("suscl1" + row + ", " + col + ", " + GlobalMethods.t);

                            for (int layupdate = 0; layupdate < GlobalMethods.max_soil_layers; layupdate++)
                            {
                                GlobalMethods.layerthickness_m[row, col, layupdate] = thickness_calc(row, col, layupdate);
                            }
                            //Debug.WriteLine("Soil after splitting");
                            //if (row == 0 & col == 0) { displaysoil(row, col); }
                            //Debug.WriteLine("suscl4");
                            //displaysoil(0, 0);

                            new_soil_mass = total_soil_mass(row, col);
                            // if (numberoflayers > GlobalMethods.max_soil_layers)
                            if (Math.Abs(old_soil_mass - new_soil_mass) > 0.00000001)
                            {
                                // this should never happen, because the data of the lowest layers have then been lost.
                                Debug.WriteLine(" Warning - loss of soil data ");
                                //displaysoil(row, col);
                            }
                            if (guiVariables.Timeseries.Timeseries_number_soil_thicker_checkbox && System.Convert.ToDouble(guiVariables.Timeseries.Timeseries_soil_thicker_textbox) < depth_m) { number_soil_thicker_than++; }
                            if (guiVariables.Timeseries.Total_average_soilthickness_checkbox) { total_average_soilthickness_m += depth_m; }
                            if (guiVariables.Timeseries.Timeseries_soil_depth_checkbox && System.Convert.ToInt32(guiVariables.Timeseries.Timeseries_soil_cell_row) == row && System.Convert.ToInt32(guiVariables.Timeseries.Timeseries_soil_cell_col) == col)
                            {
                                local_soil_depth_m = depth_m;
                            }

                            //Debug.WriteLine("suscl0" + row + ", " + col + ", " + GlobalMethods.t + " " + total_soil_mass(row, col));
                            // update GlobalMethods.dtm and soil thickness map //MMS
                            double old_thickness = GlobalMethods.soildepth_m[row, col];
                            double new_thickness = total_soil_thickness(row, col);
                            GlobalMethods.dtm[row, col] += new_thickness - old_thickness;
                            GlobalMethods.soildepth_m[row, col] = new_thickness;
                            GlobalMethods.dtmchange[row, col] += new_thickness - old_thickness;
                            GlobalMethods.dz_soil[row, col] += new_thickness - old_thickness;
                        }// end GlobalMethods.dtm!=-9999

                    } // end col
                } // end row

                if (guiVariables.Timeseries.Total_average_soilthickness_checkbox) { total_average_soilthickness_m /= GlobalMethods.number_of_data_cells; }
            }
            catch { }
            double mass_after = total_catchment_mass();
            //if (Math.Round(mass_before, 3) != Math.Round(mass_after, 3)) { Debugger.Break(); }


        }

        void soil_update_split_and_combine_layers()
        {
            double dz_standard = 0.1;
            double tolerance = 0.55; // fraction of standard thickness
                                     //where at the end of soil development, splitting and combining of soil layers is performed.

            //per layer: if too thin: combine with one of the two neighbours (the closest one in properties). 
            // too thick: split
            // if total too much - combine the most similar two layers although the product conforms to most restrictive rule about thicknesses
            //maat voor verschil is som (absolute differences in de vijf text classes and two organic matter classes)

            // 0-50 cm    min 2.5   insteek 5    maximum 10 cm       n=10    bovenste laag geen minimum (sediment HOEFT niet meteen weggemiddeld te worden - pas als nodig)
            // 50-200 cm  min 10    insteek 15    maximum 50 cm      n=8
            // daarna     min 40    insteek 50  geen max            n=4

            //for combining : only when matrix size exceeded, then those two layers that are most equal are combined
            //do not combine two layers that together are too thick for their position in the profile. If needed, make lowest layer thicker. Therefore, approach from above

            //MISSING: HOW DO WE GO FROM MAX TO LESS LAYERS? ART
            // Debug.WriteLine("suscl1 ");
            if (NA_in_map(GlobalMethods.dtm) > 0 | NA_in_map(GlobalMethods.soildepth_m) > 0)
            {
                Debug.WriteLine("err_uscl7");
            }

            //displaysoil(0,0);
            double mass_before = total_catchment_mass();
            total_average_soilthickness_m = 0;
            number_soil_thicker_than = 0;
            number_soil_coarser_than = 0;
            local_soil_depth_m = 0;
            local_soil_mass_kg = 0;

            int layer;
            int numberoflayers = 0;
            double depth_m;  // keep better track of this, currently not OK yet
            try
            {
                /*
                //displaysoil(0, 0);
                for (row = 0; row < GlobalMethods.nr; row++)
                {
                    for (col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (GlobalMethods.dtm[row, col] != -9999)
                        {
                            update_all_soil_thicknesses(row, col);
                            update_all_soil_thicknesses(row, col);
                            // Debug.WriteLine("suscl2");
                            if (NA_in_soil(row, col)) { Debugger.Break(); }


                            depth_m = 0;
                            for (layer = 0; layer < GlobalMethods.max_soil_layers; layer++)
                            {
                                ////update the layers' thickness now that textures and organic matter amounts have changed (if there is anything in the layer at all).
                                //if (!(GlobalMethods.texture_kg[row, col, layer, 0] == 0 && GlobalMethods.texture_kg[row, col, layer, 1] == 0 && GlobalMethods.texture_kg[row, col, layer, 2] == 0 && GlobalMethods.texture_kg[row, col, layer, 3] == 0 && GlobalMethods.texture_kg[row, col, layer, 4] == 0 && GlobalMethods.young_SOM_kg[row, col, layer] == 0 && GlobalMethods.old_SOM_kg[row, col, layer] == 0))
                                //{
                                //GlobalMethods.layerthickness_m[row, col, layer] = thickness_calc(row, col, layer);
                                //GlobalMethods.layerthickness_m[row, col, layer] = thickness_calc(row, col, layer);
                                //if (GlobalMethods.layerthickness_m[row,col,layer] < 0) { Debugger.Break(); } //MMS
                                //find_negative_texture_rcl(row, col, layer); //MMS
                                //}
                                if (timeseries.timeseries_soil_mass_checkbox.Checked && System.Convert.ToInt32(timeseries.timeseries_soil_cell_row.Text) == row && System.Convert.ToInt32(timeseries.timeseries_soil_cell_col.Text) == col)
                                {
                                    local_soil_mass_kg += GlobalMethods.texture_kg[row, col, layer, 0] + GlobalMethods.texture_kg[row, col, layer, 1] + GlobalMethods.texture_kg[row, col, layer, 2] + GlobalMethods.texture_kg[row, col, layer, 3] + GlobalMethods.texture_kg[row, col, layer, 4] + GlobalMethods.young_SOM_kg[row, col, layer] + GlobalMethods.old_SOM_kg[row, col, layer];
                                    if (local_soil_mass_kg < 0) { Debugger.Break(); } //MMS
                                }
                                if (timeseries.timeseries_soil_mass_checkbox.Checked && layer == 0 && GlobalMethods.texture_kg[row, col, layer, 0] / (GlobalMethods.texture_kg[row, col, layer, 0] + GlobalMethods.texture_kg[row, col, layer, 1] + GlobalMethods.texture_kg[row, col, layer, 2] + GlobalMethods.texture_kg[row, col, layer, 3] + GlobalMethods.texture_kg[row, col, layer, 4] + GlobalMethods.young_SOM_kg[row, col, layer] + GlobalMethods.old_SOM_kg[row, col, layer]) > System.Convert.ToDouble(timeseries.timeseries_soil_coarser_fraction_textbox.Text))
                                {
                                    number_soil_coarser_than++;
                                }
                            }

                        }
                    }
                } */



                //displaysoil(0, 0);
                depth_m = 0;
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (GlobalMethods.dtm[row, col] != -9999)
                        {
                            remove_empty_layers(row, col);
                            remove_empty_layers(row, col);

                            update_all_soil_thicknesses(row, col);

                            //if (NA_in_soil(row, col))
                            //{
                            //    Debug.WriteLine("err_uscl8");
                            //}


                            double old_soil_mass = total_soil_mass(row, col), new_soil_mass;
                            // Debug.WriteLine("suscl0" + row + ", " + col + ", " + GlobalMethods.t + " " + total_soil_mass(row, col));
                            //Debug.WriteLine("soil before splitting");
                            // if (row == 0 & col == 0) { displaysoil(row, col); }
                            depth_m = 0; numberoflayers = 0;
                            bool boolsplit = false;
                            bool boolcombine = false;
                            for (layer = 0; layer < (GlobalMethods.max_soil_layers - 1); layer++)
                            {
                                if (total_layer_mass(row, col, layer) > 0)
                                {

                                    //Debug.WriteLine("depth is now " + depth + " for lyr " +  layer);
                                    numberoflayers++;
                                    if (layer == 0)
                                    {
                                        if (GlobalMethods.layerthickness_m[row, col, layer] < 0.001 | total_soil_mass(row, col) < 0.001) // smaller than one mm, lighter than 1 gram -> merge with layer below, to avoid numerical problems when always a fraction leaves the profile (e.g. with GlobalMethods.creep)
                                        {
                                            combine_layers(row, col, layer, layer + 1);
                                            update_all_soil_thicknesses(row, col);
                                            boolcombine = true;
                                            if (Math.Round(old_soil_mass, 6) != Math.Round(total_soil_mass(row, col), 6))
                                            {
                                                Debug.WriteLine("err_uscl9");
                                            }
                                        }
                                        while (GlobalMethods.layerthickness_m[row, col, layer] > dz_standard * (1 + tolerance)) //Higher end, split
                                        { //split 
                                          // Debug.WriteLine("splitting after combining 1");
                                          // Debug.WriteLine("d_layer {0}", GlobalMethods.layerthickness_m[row, col, layer]);

                                            split_layer(row, col, layer, depth_m);
                                            // Debug.WriteLine("d_layer {0}", GlobalMethods.layerthickness_m[row, col, layer]);
                                            update_all_soil_thicknesses(row, col);
                                            // Debug.WriteLine("d_layer {0}", GlobalMethods.layerthickness_m[row, col, layer] );
                                            boolsplit = true;
                                        }
                                    }
                                    if (layer != 0)
                                    {
                                        if (GlobalMethods.layerthickness_m[row, col, layer] < (dz_standard * (1 - tolerance))) // Lower end, combine
                                        {
                                            if (layer_difference(row, col, layer, layer - 1) > layer_difference(row, col, layer, layer + 1))
                                            {
                                                if (Math.Abs(old_soil_mass - total_soil_mass(row, col)) > 0.00001)
                                                {
                                                    Debug.WriteLine("err_uscl10");
                                                }
                                                combine_layers(row, col, layer, layer + 1);
                                                update_all_soil_thicknesses(row, col);
                                                boolcombine = true;
                                                if (Math.Round(old_soil_mass, 6) != Math.Round(total_soil_mass(row, col), 6))
                                                {
                                                    Debug.WriteLine("err_uscl11");
                                                }
                                            }
                                            else
                                            {
                                                if (Math.Abs(old_soil_mass - total_soil_mass(row, col)) > 0.00001)
                                                {
                                                    Debug.WriteLine("err_uscl12");
                                                }
                                                combine_layers(row, col, layer - 1, layer);
                                                layer--;  //because we combined with the previous one, the current one has been replaced with one that has not yet been considered
                                                update_all_soil_thicknesses(row, col);
                                                boolcombine = true;
                                                // if (Math.Round(old_soil_mass, 6) != Math.Round(total_soil_mass(row, col), 6)) { Debugger.Break(); }
                                                if (Math.Abs(old_soil_mass - total_soil_mass(row, col)) > 0.00001)
                                                {
                                                    Debug.WriteLine("err_uscl13");
                                                }

                                            }
                                        }
                                        if (Math.Abs(old_soil_mass - total_soil_mass(row, col)) > 0.00001)
                                        {
                                            Debug.WriteLine("err_uscl14");
                                        }
                                        if (NA_in_soil(row, col))
                                        {
                                            Debug.WriteLine("err_uscl15");
                                        }

                                        while (GlobalMethods.layerthickness_m[row, col, layer] > dz_standard * (1 + tolerance)) //Higher end, split
                                        { //split 
                                          // Debug.WriteLine("splitting after combining 1");
                                            split_layer(row, col, layer, depth_m);
                                            update_all_soil_thicknesses(row, col);
                                            boolsplit = true;
                                            // if (Math.Abs(old_soil_mass-total_soil_mass(row, col))>0.00000001) { Debugger.Break(); }
                                        }
                                        if (Math.Abs(old_soil_mass - total_soil_mass(row, col)) > 0.00001)
                                        {
                                            Debug.WriteLine("err_uscl16");
                                        }

                                        //Debug.WriteLine("depth is now " + depth + " and number of layers is  " + numberoflayers);
                                    }
                                    if (Math.Abs(old_soil_mass - total_soil_mass(row, col)) > 0.00001)
                                    {
                                        Debug.WriteLine("err_uscl17");
                                    }

                                }
                            } // end layer
                              //Debug.WriteLine("suscl1" + row + ", " + col + ", " + GlobalMethods.t);
                            if (NA_in_soil(row, col))
                            {
                                Debug.WriteLine("err_uscl18");
                            }

                            for (int layupdate = 0; layupdate < GlobalMethods.max_soil_layers; layupdate++)
                            {
                                GlobalMethods.layerthickness_m[row, col, layupdate] = thickness_calc(row, col, layupdate);
                            }
                            //Debug.WriteLine("Soil after splitting");
                            //if (row == 0 & col == 0) { displaysoil(row, col); }

                            // Debug.WriteLine("suscl4");
                            //displaysoil(0, 0);
                            if (Math.Abs(old_soil_mass - total_soil_mass(row, col)) > 0.00001)
                            {
                                Debug.WriteLine("err_uscl19");
                            }
                            new_soil_mass = total_soil_mass(row, col);
                            // if (numberoflayers > GlobalMethods.max_soil_layers)
                            if (Math.Abs(old_soil_mass - new_soil_mass) > 0.00000001)
                            {
                                // this should never happen, because the data of the lowest layers have then been lost.
                                // Debug.WriteLine("{0}", GlobalMethods.t);
                                Debug.WriteLine(" Warning - loss of soil data ");
                                //displaysoil(row, col);
                                Debug.WriteLine("err_uscl20");

                            }
                            if (guiVariables.Timeseries.Timeseries_number_soil_thicker_checkbox && System.Convert.ToDouble(guiVariables.Timeseries.Timeseries_soil_thicker_textbox) < depth_m) { number_soil_thicker_than++; }
                            if (guiVariables.Timeseries.Total_average_soilthickness_checkbox) { total_average_soilthickness_m += depth_m; }
                            if (guiVariables.Timeseries.Timeseries_soil_depth_checkbox && System.Convert.ToInt32(guiVariables.Timeseries.Timeseries_soil_cell_row) == row && System.Convert.ToInt32(guiVariables.Timeseries.Timeseries_soil_cell_col) == col)
                            {
                                local_soil_depth_m = depth_m;
                            }

                            //Debug.WriteLine("suscl0" + row + ", " + col + ", " + GlobalMethods.t + " " + total_soil_mass(row, col));
                            // update GlobalMethods.dtm and soil thickness map //MMS
                            update_all_soil_thicknesses(row, col);
                            double old_thickness = GlobalMethods.soildepth_m[row, col];
                            double new_thickness = total_soil_thickness(row, col);
                            GlobalMethods.dtm[row, col] += new_thickness - old_thickness;
                            GlobalMethods.soildepth_m[row, col] = new_thickness;
                            GlobalMethods.dtmchange[row, col] += new_thickness - old_thickness;
                            GlobalMethods.dz_soil[row, col] += new_thickness - old_thickness;
                        } // end GlobalMethods.dtm!=-9999
                    } // end col
                } // end row

                if (guiVariables.Timeseries.Total_average_soilthickness_checkbox) { total_average_soilthickness_m /= GlobalMethods.number_of_data_cells; }
            }
            catch
            {
                Debug.WriteLine("err_uscl21");
            }
            double mass_after = total_catchment_mass();

        } // aangepast voor constante diktes


        private void find_negative_texture()
        {
            for (int rr = 0; rr < GlobalMethods.nr; rr++)
            {
                for (int cc = 0; cc < GlobalMethods.nc; cc++)
                {
                    for (int ll = 0; ll < GlobalMethods.max_soil_layers; ll++)
                    {
                        for (int tt = 0; tt < 5; tt++)
                        {
                            if (GlobalMethods.texture_kg[rr, cc, ll, tt] < 0)
                            {
                                Debug.WriteLine("err_nt1");
                            }
                        }
                    }
                }
            }
        }

        double layer_difference(int rowwer, int coller, int lay1, int lay2)   //calculates a simple measure of difference between two soil layers based on the sum of relative differences in a set of properties
        {
            double average_property_value = 0, property_difference = 0, sum_property_difference = 0;

            try
            {
                // double average_property_value = 0, property_difference = 0, sum_property_difference = 0;
                for (i = 0; i < 5; i++)
                {
                    //account for total soil mass
                    average_property_value = (GlobalMethods.texture_kg[rowwer, coller, lay1, i] + GlobalMethods.texture_kg[rowwer, coller, lay2, i]) / 2;
                    property_difference = Math.Abs(GlobalMethods.texture_kg[rowwer, coller, lay1, i] - GlobalMethods.texture_kg[rowwer, coller, lay2, i]);
                    sum_property_difference += property_difference / average_property_value;
                }
                average_property_value = (GlobalMethods.old_SOM_kg[rowwer, coller, lay1] + GlobalMethods.old_SOM_kg[rowwer, coller, lay2]) / 2;
                property_difference = Math.Abs(GlobalMethods.old_SOM_kg[rowwer, coller, lay1] - GlobalMethods.old_SOM_kg[rowwer, coller, lay2]);
                sum_property_difference += property_difference / average_property_value;
                average_property_value = (GlobalMethods.young_SOM_kg[rowwer, coller, lay1] + GlobalMethods.young_SOM_kg[rowwer, coller, lay2]) / 2;
                property_difference = Math.Abs(GlobalMethods.young_SOM_kg[rowwer, coller, lay1] - GlobalMethods.young_SOM_kg[rowwer, coller, lay2]);
                sum_property_difference += property_difference / average_property_value;
                sum_property_difference /= 7;
                return sum_property_difference;
            }
            catch
            {
                MessageBox.Show("error in calculating layer difference");
            }
            return sum_property_difference;
        }

        void combine_layers(int rowwer, int coller, int lay1, int lay2)  // combines two soil layers into the first, recalculates their new thickness and shifts underlying layers up
        {

            double mass_before = total_layer_mass(rowwer, coller, lay1) + total_layer_mass(rowwer, coller, lay2);
            try
            {
                double old_soil_mass1 = total_soil_mass(rowwer, coller);
                // Debug.WriteLine("Total soil mass: {0}", old_soil_mass); displaysoil(rowwer, coller); 
                // Debug.WriteLine("cl0");
                //Debug.WriteLine("total soil mass = " + total_soil_mass(rowwer, coller));
                for (i = 0; i < 5; i++)
                {
                    GlobalMethods.texture_kg[rowwer, coller, lay1, i] += GlobalMethods.texture_kg[rowwer, coller, lay2, i];
                    GlobalMethods.texture_kg[rowwer, coller, lay2, i] = 0;// set to zero. otherwise the shifting of the layers doesn'GlobalMethods.t work
                    double new_soil_mass = total_soil_mass(rowwer, coller);

                }

                GlobalMethods.old_SOM_kg[rowwer, coller, lay1] += GlobalMethods.old_SOM_kg[rowwer, coller, lay2];
                GlobalMethods.old_SOM_kg[rowwer, coller, lay2] = 0;

                GlobalMethods.young_SOM_kg[rowwer, coller, lay1] += GlobalMethods.young_SOM_kg[rowwer, coller, lay2];
                GlobalMethods.young_SOM_kg[rowwer, coller, lay2] = 0;

                //Debug.WriteLine("cl1");
                //Debug.WriteLine("total soil mass = " + total_soil_mass(rowwer, coller));
                GlobalMethods.layerthickness_m[rowwer, coller, lay1] = thickness_calc(rowwer, coller, lay1);    // thickness_calc uses a pdf to calculate bulk density and hence layer thickness

                for (int layert = lay2; layert < GlobalMethods.max_soil_layers - 1; layert++) // for all underlying layers, shift one up (since there is space anyway)
                {
                    // Debug.WriteLine(layert);
                    for (i = 0; i < 5; i++)
                    {
                        GlobalMethods.texture_kg[rowwer, coller, layert, i] = GlobalMethods.texture_kg[rowwer, coller, layert + 1, i];
                    }
                    GlobalMethods.old_SOM_kg[rowwer, coller, layert] = GlobalMethods.old_SOM_kg[rowwer, coller, layert + 1];
                    GlobalMethods.young_SOM_kg[rowwer, coller, layert] = GlobalMethods.young_SOM_kg[rowwer, coller, layert + 1];
                }

                //now set the last layer to sentinel value of -1
                for (i = 0; i < 5; i++)
                {
                    GlobalMethods.texture_kg[rowwer, coller, GlobalMethods.max_soil_layers - 1, i] = 0;
                }
                GlobalMethods.old_SOM_kg[rowwer, coller, GlobalMethods.max_soil_layers - 1] = 0;
                GlobalMethods.young_SOM_kg[rowwer, coller, GlobalMethods.max_soil_layers - 1] = 0;
                GlobalMethods.layerthickness_m[rowwer, coller, GlobalMethods.max_soil_layers - 1] = 0;
                GlobalMethods.bulkdensity[rowwer, coller, GlobalMethods.max_soil_layers - 1] = 0;

                double new_soil_mass1 = total_soil_mass(rowwer, coller);
                //if (Math.Abs(old_soil_mass1 - new_soil_mass1) > 0.00000001) { displaysoil(rowwer, coller); Debugger.Break(); }
                double mass_after = total_layer_mass(rowwer, coller, lay1);
                if (Math.Abs(mass_before - mass_after) > 0.0001)
                {
                    Debug.WriteLine("err_cl1");
                }

            }
            catch
            {
                Debug.WriteLine("Failed at combining layer at row {0}, col {1} at time {2}", row, col, GlobalMethods.t);
            }
        }

        public double bulk_density_calc(double coarse_mass, double sand_mass, double silt_mass, double clay_mass, double fine_clay_mass, double OMo_mass, double OMy_mass, double depth)
        {
            if (depth == 0) { depth = 0.001; } //MvdM value to prevent no data when calculating BD for layers that were previously empty, or for sediments that are deposited. A value of 1 results in an infinite value
            double bd = 2700, combined_frac, m_finesoil;
            m_finesoil = sand_mass + silt_mass + clay_mass + fine_clay_mass;
            if (m_finesoil > 0)
            {
                combined_frac = sand_mass / m_finesoil + 0.76 * silt_mass / m_finesoil;

                bd = 1000 * (1.35 + 0.00452 * 100 * combined_frac + Math.Pow(100 * combined_frac - 44.65, 2) * -0.0000614 + 0.06 * Math.Log10(depth));  // in kg/m3

                //now coarse fragment and SOM correction
                bd = (coarse_mass + m_finesoil + OMo_mass + OMy_mass) / ((m_finesoil / bd) + (coarse_mass / 2700) + (OMo_mass + OMy_mass) / 224); // ooit through interface    
            }
            else
            {
                bd = 2700;
            }
            if (double.IsNaN(bd) | double.IsInfinity(bd)) { Debugger.Break(); }
            return bd;
            // return 1500;
        }

        double thickness_calc(int rowwer, int coller, int lay1)
        {
            double thickness, soil_mass = 0;

            // calculate current depth of layer, for bulk density calculations, using current thickness. 
            double depth_m = 0;
            for (int lay_temp = 0; lay_temp < lay1; lay_temp++)
            {
                depth_m += GlobalMethods.layerthickness_m[rowwer, coller, lay_temp];
            }
            depth_m += GlobalMethods.layerthickness_m[rowwer, coller, lay1] / 2;

            int i;
            //first calculate total soil mass to calculate mass percentages for the size fractions
            for (i = 1; i < 5; i++)
            {
                soil_mass += GlobalMethods.texture_kg[rowwer, coller, lay1, i];
            }
            soil_mass += GlobalMethods.old_SOM_kg[rowwer, coller, lay1] + GlobalMethods.young_SOM_kg[rowwer, coller, lay1];
            //if (GlobalMethods.t == 6 && row == 193 && col == 58) { Debug.WriteLine("A " + row + " " + col + " soil mass " + soil_mass); displaysoil(row, col); }
            if (soil_mass > 0)
            {
                GlobalMethods.bulkdensity[rowwer, coller, lay1] = bulk_density_calc(GlobalMethods.texture_kg[rowwer, coller, lay1, 0], GlobalMethods.texture_kg[rowwer, coller, lay1, 1], GlobalMethods.texture_kg[rowwer, coller, lay1, 2], GlobalMethods.texture_kg[rowwer, coller, lay1, 3], GlobalMethods.texture_kg[rowwer, coller, lay1, 4], GlobalMethods.old_SOM_kg[rowwer, coller, lay1], GlobalMethods.young_SOM_kg[rowwer, coller, lay1], depth_m);
                //sand_fraction = GlobalMethods.texture_kg[rowwer, coller, lay1, 1] / soil_mass;
                //silt_fraction = GlobalMethods.texture_kg[rowwer, coller, lay1, 2] / soil_mass;
                //combined_fraction = sand_fraction + 0.76 * silt_fraction;
                //GlobalMethods.bulkdensity[rowwer, coller, lay1] = 1000 * (1.35 + 0.00452 * 100 * combined_fraction + Math.Pow(100 * combined_fraction - 44.65, 2) * -0.0000614);  // in kg/m3
                ////now coarse fragment correction

                //GlobalMethods.bulkdensity[rowwer, coller, lay1] = (soil_mass + GlobalMethods.texture_kg[rowwer, coller, lay1, 0]) / ((soil_mass / GlobalMethods.bulkdensity[rowwer, coller, lay1]) + (GlobalMethods.texture_kg[rowwer, coller, lay1, 0] / 2700)); // ooit through interface              
                // if (GlobalMethods.t == 6 && row == 193 && col == 58) { Debug.WriteLine("B " + row + " " + col + " soil mass " + soil_mass); }
            }
            else
            {
                //either there is no soil at all, or there is only coarse material
                if (GlobalMethods.texture_kg[rowwer, coller, lay1, 0] > 0)
                {
                    GlobalMethods.bulkdensity[rowwer, coller, lay1] = 2700;   //kg/m3
                }
            }
            // if (GlobalMethods.t == 6 && row == 193 && col == 58) { Debug.WriteLine("C " + row + " " + col + " soil mass " + soil_mass); }
            soil_mass += GlobalMethods.texture_kg[rowwer, coller, lay1, 0];
            if (soil_mass == 0)
            {
                thickness = 0;
            }
            else
            {
                thickness = (soil_mass) / (GlobalMethods.dx * GlobalMethods.dx * GlobalMethods.bulkdensity[rowwer, coller, lay1]);  // thickness in m per unit area
            }
            if (double.IsNaN(thickness)) { thickness = 0.00000000001; }
            return thickness;
        }

        double total_soil_mass(int rowmass, int colmass)
        {
            double tot_mass = 0;
            for (int lay = 0; lay < GlobalMethods.max_soil_layers; lay++)
            {
                for (int ii = 0; ii < GlobalMethods.n_texture_classes; ii++)
                {
                    tot_mass += GlobalMethods.texture_kg[rowmass, colmass, lay, ii];
                }
                tot_mass += GlobalMethods.old_SOM_kg[rowmass, colmass, lay];
                tot_mass += GlobalMethods.young_SOM_kg[rowmass, colmass, lay];
            }
            return (tot_mass);
        }

        double total_layer_mass(int rowmass, int colmass, int laymass)
        {
            double tot_mass = 0;

            for (int ii = 0; ii < 5; ii++)
            {
                tot_mass += GlobalMethods.texture_kg[rowmass, colmass, laymass, ii];
            }
            tot_mass += GlobalMethods.old_SOM_kg[rowmass, colmass, laymass];
            tot_mass += GlobalMethods.young_SOM_kg[rowmass, colmass, laymass];

            return (tot_mass);
        }

        double total_layer_fine_earth_mass(int rowmass, int colmass, int laymass)
        {
            double tot_mass = 0;

            for (int ii = 1; ii < 5; ii++)
            {
                tot_mass += GlobalMethods.texture_kg[rowmass, colmass, laymass, ii];
            }
            tot_mass += GlobalMethods.old_SOM_kg[rowmass, colmass, laymass];
            tot_mass += GlobalMethods.young_SOM_kg[rowmass, colmass, laymass];

            return (tot_mass);
        }

        bool findnegativetexture()
        {
            bool neg = false;

            try
            {
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        for (int lay = 0; lay < GlobalMethods.max_soil_layers; lay++)
                        {
                            for (int tex = 0; tex < GlobalMethods.n_texture_classes; tex++)
                            {
                                if (GlobalMethods.texture_kg[row, col, lay, tex] < 0) { neg = true; }
                            }
                        }
                    }
                }
            }
            catch
            {
                Debug.WriteLine("err_nt2");

            }

            return neg;
        }

        double total_mass_in_transport()
        {
            double tot_mass = 0;
            for (int rowmass = 0; rowmass < GlobalMethods.nr; rowmass++)
            {
                for (int colmass = 0; colmass < GlobalMethods.nc; colmass++)
                {
                    for (int ii = 0; ii < 5; ii++)
                    {
                        tot_mass += GlobalMethods.sediment_in_transport_kg[rowmass, colmass, ii];
                    }

                    tot_mass += GlobalMethods.old_SOM_in_transport_kg[rowmass, colmass];
                    tot_mass += GlobalMethods.young_SOM_in_transport_kg[rowmass, colmass];
                }
            }
            return (tot_mass);
        }

        double mass_in_transport_row_col(int row1, int col1)
        {
            double tot_mass = 0;
            for (int ii = 0; ii < 5; ii++)
            {
                tot_mass += GlobalMethods.sediment_in_transport_kg[row1, col1, ii];
            }

            tot_mass += GlobalMethods.old_SOM_in_transport_kg[row1, col1];
            tot_mass += GlobalMethods.young_SOM_in_transport_kg[row1, col1];

            return (tot_mass);
        }

        double total_catchment_mass()
        {
            double tot_mass = 0;
            for (int rowmass = 0; rowmass < GlobalMethods.nr; rowmass++)
            {
                for (int colmass = 0; colmass < GlobalMethods.nc; colmass++)
                {
                    for (int lay = 0; lay < GlobalMethods.max_soil_layers; lay++)
                    {
                        for (int ii = 0; ii < 5; ii++)
                        {
                            tot_mass += GlobalMethods.texture_kg[rowmass, colmass, lay, ii];
                        }
                        tot_mass += GlobalMethods.old_SOM_kg[rowmass, colmass, lay];
                        tot_mass += GlobalMethods.young_SOM_kg[rowmass, colmass, lay];
                    }
                }
            }


            return (tot_mass);
        }

        double total_catchment_elevation()
        {
            double tot_elev = 0;
            for (int rowmass = 0; rowmass < GlobalMethods.nr; rowmass++)
            {
                for (int colmass = 0; colmass < GlobalMethods.nc; colmass++)
                {
                    if (GlobalMethods.dtm[rowmass, colmass] != -9999)
                    {
                        tot_elev += GlobalMethods.dtm[rowmass, colmass];
                    }

                }
            }


            return (tot_elev);
        }

        double total_soil_thickness(int rowthick, int colthick)
        {
            double tot_thick = 0;
            for (int lay = 0; lay < GlobalMethods.max_soil_layers; lay++)
            {
                tot_thick += GlobalMethods.layerthickness_m[rowthick, colthick, lay];
            }
            return (tot_thick);
        }

        void split_layer(int rowwer, int coller, int lay1, double currentdepth) // splits layers 
        {
            try
            {
                double max_layer_difference, current_difference, maximum_allowed_thickness;
                //splitting will increase the number of layers. If this splits beyond the max number of layers, then combine the two most similar ones 
                int laynum, combininglayer = -1;
                // Debug.WriteLine("sl0");
                //Debug.WriteLine("total soil mass = " + total_soil_mass(rowwer, coller));
                double mass_lowest_layer = total_layer_mass(rowwer, coller, GlobalMethods.max_soil_layers - 1);
                // Debug.WriteLine("total mass last layer = {0}", mass_lowest_layer);
                if ((total_layer_mass(rowwer, coller, GlobalMethods.max_soil_layers - 1) > 0))  // so, if we are using the lowest possible layer already:
                {
                    //if they are already all in use, then the split will create one too many. We start by looking for the two most similar layers that would not create a too-thick product (do we need to do that last part?)
                    max_layer_difference = 100; //100 is a huge difference
                    for (laynum = 0; laynum < GlobalMethods.max_soil_layers - 1; laynum++)
                    {
                        current_difference = layer_difference(rowwer, coller, laynum, laynum + 1);
                        maximum_allowed_thickness = 9999;   // 9999 is a sentinel value and means infinitely thick 
                        if (currentdepth < 0.5) { maximum_allowed_thickness = 0.1; }
                        else { if (currentdepth < 2.5) { maximum_allowed_thickness = 0.5; } }
                        maximum_allowed_thickness = dz_standard * (1 + tolerance);
                        if (GlobalMethods.layerthickness_m[rowwer, coller, laynum] + GlobalMethods.layerthickness_m[rowwer, coller, laynum + 1] < maximum_allowed_thickness)  //if it potentially is possible to combine them 
                        {
                            if (current_difference <= max_layer_difference)   // the equal to condition means that we prefer to combine layers lower in the profile (if equally different from each other)
                            {
                                max_layer_difference = current_difference; combininglayer = laynum;
                            }
                        }
                    }
                    //Debug.WriteLine("sl1");
                    //Debug.WriteLine("total soil mass = " + total_soil_mass(rowwer, coller));

                    //combine the two most-similar layers, or the lowest two if nothing else possible
                    if (max_layer_difference == 100) { combininglayer = GlobalMethods.max_soil_layers - 2; }
                    //Debug.WriteLine("test");
                    //Debug.WriteLine(rowwer + "," + coller + "," + combininglayer + "," + combininglayer + 1+","+GlobalMethods.t);
                    try { combine_layers(rowwer, coller, combininglayer, combininglayer + 1); }
                    catch { Debug.WriteLine(" failed to combine layers to prepare for splitting "); }
                    //make sure to change lay1 if needed (because something overlying has been combined into 1 for instance).
                    if (combininglayer == lay1 || combininglayer == lay1 - 1)
                    {
                        Debug.WriteLine("the layer that needed to be split has now been combined: layer {0} at GlobalMethods.t {1}", combininglayer, GlobalMethods.t);
                        // displaysoil(rowwer, coller);
                        // Debugger.Break();
                    }
                    if (combininglayer < lay1) { lay1++; } // mvdm changed from++ to --, because if layers above are combined, the target layer has moved up one spot
                }
                //Debug.WriteLine("sl2");
                //Debug.WriteLine("total soil mass = " + total_soil_mass(rowwer, coller));

                // now we can move all layers down below the one we want to split
                for (laynum = GlobalMethods.max_soil_layers - 1; laynum >= lay1 + 2; laynum--)  // we want to clear layer lay1+1 (so we run through move-receiving layers from below up to lay1+2). 
                                                                                  //This means that layer laynum+1, into which we want to split, will be evacuated and will give its values to laynum+2;
                {
                    // Debug.WriteLine("sl2a, laynum = "+laynum+", lay1+2 = "+lay1 + 2+"tex lay 19 =" + GlobalMethods.texture_kg[rowwer,coller,19,2]);
                    for (i = 0; i < 5; i++)
                    {
                        GlobalMethods.texture_kg[rowwer, coller, laynum, i] = GlobalMethods.texture_kg[rowwer, coller, laynum - 1, i];
                    }
                    // Debug.WriteLine("sl2b, ");
                    GlobalMethods.old_SOM_kg[rowwer, coller, laynum] = GlobalMethods.old_SOM_kg[rowwer, coller, laynum - 1];
                    GlobalMethods.young_SOM_kg[rowwer, coller, laynum] = GlobalMethods.young_SOM_kg[rowwer, coller, laynum - 1];
                    // OSL_age[rowwer, coller, laynum] = OSL_age[rowwer, coller, laynum - 1];
                }
                //Debug.WriteLine("sl3");
                //Debug.WriteLine("total soil mass = " + total_soil_mass(rowwer, coller));

                //Debug.WriteLine(" moved layers one down to make space for split layer ");
                if ((lay1 + 1) == (GlobalMethods.max_soil_layers - 1))
                {
                    double div = 0.1 / (GlobalMethods.layerthickness_m[rowwer, coller, lay1]); // aim to have the split layer at 0.1 m
                    if (div > 1) { div = 1; }
                    for (i = 0; i < 5; i++)
                    {
                        GlobalMethods.texture_kg[rowwer, coller, lay1 + 1, i] += GlobalMethods.texture_kg[rowwer, coller, lay1, i] * (1 - div); GlobalMethods.texture_kg[rowwer, coller, lay1, i] *= div;
                        if (double.IsNaN(GlobalMethods.texture_kg[rowwer, coller, lay1, i]))
                        {
                            Debug.WriteLine("err_spl1");
                        }
                        if (double.IsNaN(GlobalMethods.texture_kg[rowwer, coller, lay1 + 1, i]))
                        {
                            Debug.WriteLine("err_spl2");
                        }
                    }
                    GlobalMethods.old_SOM_kg[rowwer, coller, lay1 + 1] += GlobalMethods.old_SOM_kg[rowwer, coller, lay1] * (1 - div); GlobalMethods.old_SOM_kg[rowwer, coller, lay1] *= div;
                    GlobalMethods.young_SOM_kg[rowwer, coller, lay1 + 1] += GlobalMethods.young_SOM_kg[rowwer, coller, lay1] * (1 - div); GlobalMethods.young_SOM_kg[rowwer, coller, lay1] *= div;
                    //Debug.WriteLine(" successfully split layer ");
                    //Debug.WriteLine("sl4");
                    //Debug.WriteLine("total soil mass = " + total_soil_mass(rowwer, coller));
                }
                else // even splitting
                {
                    for (i = 0; i < 5; i++)
                    {
                        GlobalMethods.texture_kg[rowwer, coller, lay1 + 1, i] = GlobalMethods.texture_kg[rowwer, coller, lay1, i] / 2; GlobalMethods.texture_kg[rowwer, coller, lay1, i] /= 2;
                    }
                    GlobalMethods.old_SOM_kg[rowwer, coller, lay1 + 1] = GlobalMethods.old_SOM_kg[rowwer, coller, lay1] / 2; GlobalMethods.old_SOM_kg[rowwer, coller, lay1] /= 2;
                    GlobalMethods.young_SOM_kg[rowwer, coller, lay1 + 1] = GlobalMethods.young_SOM_kg[rowwer, coller, lay1] / 2; GlobalMethods.young_SOM_kg[rowwer, coller, lay1] /= 2;
                    //Debug.WriteLine(" successfully split layer ");
                    //Debug.WriteLine("sl4");
                    //Debug.WriteLine("total soil mass = " + total_soil_mass(rowwer, coller));
                }

            }
            catch
            {
                Debug.WriteLine("Failed at splitting layer at row {0}, col {1} at time {2}", row, col, GlobalMethods.t);
            }
        }

        void split_layer_till(int rowwer, int coller, int lay1, double currentdepth) // splits layers 
        {

            double max_layer_difference, current_difference, maximum_allowed_thickness, frac_ap, frac_soil;
            //splitting will increase the number of layers. If this splits beyond the max number of layers, then combine the two most similar ones 
            int laynum, combininglayer = -1;
            if ((GlobalMethods.layerthickness_m[rowwer, coller, GlobalMethods.max_soil_layers - 1] > 0))  // so, if we are using the lowest possible layer already:
            {
                //if they are already all in use, then the split will create one too many. We start by looking for the two most similar layers that would not create a too-thick product (do we need to do that last part?)
                max_layer_difference = 100; //100 is a huge difference
                for (laynum = 0; laynum < GlobalMethods.max_soil_layers - 1; laynum++)
                {
                    current_difference = layer_difference(rowwer, coller, laynum, laynum + 1);
                    maximum_allowed_thickness = 9999;   // 9999 is a sentinel value and means infinitely thick 
                    if (currentdepth < 0.5) { maximum_allowed_thickness = 0.1; }
                    else { if (currentdepth < 2.5) { maximum_allowed_thickness = 0.5; } }

                    if (GlobalMethods.layerthickness_m[rowwer, coller, laynum] + GlobalMethods.layerthickness_m[rowwer, coller, laynum + 1] < maximum_allowed_thickness)  //if it potentially is possible to combine them 
                    {
                        if (current_difference <= max_layer_difference)   // the equal to condition means that we prefer to combine layers lower in the profile (if equally different from each other)
                        {
                            max_layer_difference = current_difference; combininglayer = laynum;
                        }
                    }
                }
                //combine the two most-similar layers, or the lowest two if nothing else possible
                if (max_layer_difference == 100) { combininglayer = GlobalMethods.max_soil_layers - 2; }
                try { combine_layers(rowwer, coller, combininglayer, combininglayer + 1); }
                catch { Debug.WriteLine(" failed to combine layers to prepare for splitting "); }
                //make sure to change lay1 if needed (because something overlying has been combined into 1 for instance).
                if (combininglayer == lay1 || combininglayer == lay1 - 1) { Debug.WriteLine("the layer that needed to be split has now been combined "); }
                if (combininglayer < lay1) { lay1++; }
            }
            // now we can move all layers down below the one we want to split
            for (laynum = GlobalMethods.max_soil_layers - 1; laynum <= lay1 + 2; laynum--)  // we want to clear layer lay1+1 (so we run through move-receiving layers from below up to lay1+2). 
                                                                              //This means that layer laynum+1, into which we want to split, will be evacuated and will give its values to laynum+2;
            {
                for (i = 0; i < 5; i++)
                {
                    GlobalMethods.texture_kg[rowwer, coller, laynum, i] = GlobalMethods.texture_kg[rowwer, coller, laynum - 1, i];
                }
                GlobalMethods.old_SOM_kg[rowwer, coller, laynum] = GlobalMethods.old_SOM_kg[rowwer, coller, laynum - 1];
                GlobalMethods.young_SOM_kg[rowwer, coller, laynum] = GlobalMethods.young_SOM_kg[rowwer, coller, laynum - 1];
            }
            //Debug.WriteLine(" moved layers one down to make space for split layer ");
            double inithick = GlobalMethods.layerthickness_m[rowwer, coller, laynum];
            frac_ap = plough_depth / inithick;
            frac_soil = 1 - frac_ap;
            for (i = 0; i < 5; i++)
            {
                GlobalMethods.texture_kg[rowwer, coller, lay1 + 1, i] = GlobalMethods.texture_kg[rowwer, coller, lay1, i] * frac_soil; GlobalMethods.texture_kg[rowwer, coller, lay1, i] *= frac_ap;
            }
            GlobalMethods.old_SOM_kg[rowwer, coller, lay1 + 1] = GlobalMethods.old_SOM_kg[rowwer, coller, lay1] / 2; GlobalMethods.old_SOM_kg[rowwer, coller, lay1] /= 2;
            GlobalMethods.young_SOM_kg[rowwer, coller, lay1 + 1] = GlobalMethods.young_SOM_kg[rowwer, coller, lay1] / 2; GlobalMethods.young_SOM_kg[rowwer, coller, lay1] /= 2;
            //Debug.WriteLine(" successfully split layer ");
        }

        bool search_nodataneighbour(int row, int col)
        {
            bool ndn = false;
            for (i = (-1); i <= 1; i++)
            {
                for (j = (-1); j <= 1; j++)
                {
                    if ((((row + i) < 0) | ((row + i) >= GlobalMethods.nr)) | (((col + j) < 0) | ((col + j) >= GlobalMethods.nc)))
                    {
                        ndn = true;
                    }
                }
            }
            return (ndn);
        }

        void soil_physical_weathering()  //calculate physical weathering
        {
            double old_mass = total_catchment_mass();
            int cells = GlobalMethods.nr * GlobalMethods.nc;
            int layer, tex_class;
            double depth;
            try
            {
                //another parallelization opportunity
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        int tempcol = col;
                        depth = 0;
                        for (layer = 0; layer < GlobalMethods.max_soil_layers; layer++)
                        {
                            if (GlobalMethods.layerthickness_m[row, tempcol, layer] > 0)
                            {
                                int templayer = layer;
                                depth += GlobalMethods.layerthickness_m[row, tempcol, templayer] / 2;
                                for (tex_class = 0; tex_class <= 2; tex_class++)   //we only physically weather the coarse, sand and silt fractions.
                                {
                                    int tempclass = tex_class;
                                    // calculate the mass involved in physical weathering
                                    weathered_mass_kg = GlobalMethods.texture_kg[row, tempcol, templayer, tempclass] * physical_weathering_constant * Math.Exp(-Cone * depth) * -Ctwo / Math.Log10(upper_particle_size[tempclass]) * dt;
                                    total_phys_weathered_mass_kg += weathered_mass_kg;
                                    //Debug.WriteLine(" weathered mass is " + weathered_mass + " for class " + tempclass );
                                    // calculate the products involved
                                    if (tex_class == 0)
                                    {
                                        GlobalMethods.texture_kg[row, tempcol, templayer, tempclass + 1] += 0.975 * weathered_mass_kg;
                                        GlobalMethods.texture_kg[row, tempcol, templayer, tempclass + 2] += 0.025 * weathered_mass_kg;
                                    }
                                    if (tex_class == 1)
                                    {
                                        GlobalMethods.texture_kg[row, tempcol, templayer, tempclass + 1] += 0.96 * weathered_mass_kg;
                                        GlobalMethods.texture_kg[row, tempcol, templayer, tempclass + 2] += 0.04 * weathered_mass_kg;
                                    }
                                    if (tex_class == 2)
                                    {
                                        GlobalMethods.texture_kg[row, tempcol, templayer, tempclass + 1] += weathered_mass_kg;
                                    }
                                    GlobalMethods.texture_kg[row, tempcol, templayer, tempclass] -= weathered_mass_kg;
                                }
                                depth += GlobalMethods.layerthickness_m[row, tempcol, templayer] / 2;
                            }
                        }
                    }  //);
                } // end for cells
                  //timeseries
                if (guiVariables.Timeseries.Timeseries_cell_waterflow_check)
                {
                    timeseries_matrix[GlobalMethods.t, timeseries_order[15]] = total_phys_weathered_mass_kg;
                }
            }
            catch { Debug.WriteLine(" Soil physical weathering calculation threw an exception"); }

            double new_mass = total_catchment_mass();
            if (Math.Abs(old_mass - new_mass) > 0.000001)
            {
                Debug.WriteLine("err_spw1");
            }
        }

        void SPITS_soil_physical_weathering()  //calculate sedimentary rock (siltstone, limestoneF physical weathering
        {
            // in this variant, coarse material (siltstone, limestone) weathers only into a silt fraction. 
            // Nothing else weathers. 
            // Not all of the coarse material weathers into silt, a certain fraction is lost to dissolution (90%).

            int cells = GlobalMethods.nr * GlobalMethods.nc;
            int layer, tex_class;
            double depth;
            try
            {
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)

                    //Parallel.For(0, GlobalMethods.nc-1, col =>                    //we should paralellize over cols. Problem so far seems to be that the GlobalMethods.nc-1 or layer limit is exceeded
                    {
                        int tempcol = col;
                        //Code here is executed in parallel as much as possible for different soils in different places. 
                        //Main assumption: soils affect each other only through their surface interactions and not e.g. through throughflow
                        depth = 0;
                        for (layer = 0; layer < GlobalMethods.max_soil_layers; layer++)
                        {
                            if (GlobalMethods.layerthickness_m[row, tempcol, layer] > 0)
                            {
                                int templayer = layer;
                                depth += GlobalMethods.layerthickness_m[row, tempcol, templayer] / 2;
                                tex_class = 0;
                                int tempclass = tex_class;
                                // calculate the mass involved in physical weathering
                                weathered_mass_kg = GlobalMethods.texture_kg[row, tempcol, templayer, tempclass] * physical_weathering_constant * Math.Exp(-Cone * depth) * -Ctwo / Math.Log10(upper_particle_size[tempclass]) * dt;
                                //Debug.WriteLine(" weathered mass is " + weathered_mass + " for class " + tempclass );
                                // calculate the products involved
                                GlobalMethods.texture_kg[row, tempcol, templayer, tempclass + 2] += 0.1 * weathered_mass_kg;
                                GlobalMethods.texture_kg[row, tempcol, templayer, tempclass] -= weathered_mass_kg;
                            }
                        }
                    }  //);
                }
            }
            catch { Debug.WriteLine(" Soil physical weathering calculation threw an exception"); }

        }

        void SPITS_aeolian_deposition()
        {
            //tricks the deposition process by playing with tillage fields. Tillage shoudl be ON - but with zero par values.
            try
            {
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (GlobalMethods.tillfields[row, col] == 1)
                        {
                            //deposition in kg/m2/y is 0.063 
                            GlobalMethods.texture_kg[row, col, 0, 1] += 0.063 * GlobalMethods.dx * GlobalMethods.dx;
                        }
                    }
                }

            }
            catch { }
        }

        void soil_chemical_weathering()

        {
            int cells = GlobalMethods.nr * GlobalMethods.nc;
            int layer, tex_class;
            double depth, weathered_mass_kg, total_weath_mass, fraction_neoform;
            total_chem_weathered_mass_kg = 0;
            total_fine_neoformed_mass_kg = 0;
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    //Main assumption: soils affect each other only through their surface interactions and not e.g. through throughflow
                    depth = 0; total_weath_mass = 0;
                    for (layer = 0; layer < GlobalMethods.max_soil_layers; layer++)
                    {
                        if (GlobalMethods.layerthickness_m[row, col, layer] > 0)
                        {
                            depth += GlobalMethods.layerthickness_m[row, col, layer] / 2;
                            for (tex_class = 1; tex_class <= 4; tex_class++) // only sand, silt, clay and fine clay are chemically weathered
                            {
                                weathered_mass_kg = GlobalMethods.texture_kg[row, col, layer, tex_class] * chemical_weathering_constant / 10 * Math.Exp(-Cthree * depth) * Cfour * specific_area[tex_class] * dt;

                                if (guiVariables.Daily_water) { weathered_mass_kg *= waterfactor[row, col]; }

                                //Debug.WriteLine(" weath mass for layer " + layer + " class " + tex_class + " is " + weathered_mass_kg + " " + Math.Exp(-Cthree * depth));
                                // note that the chem_weath constant is in kg/m2 mineral surface / y (in contrast to the original value from Salvador Blanes (mol/m2 mineral/s)
                                if (weathered_mass_kg > GlobalMethods.texture_kg[row, col, layer, tex_class]) { weathered_mass_kg = GlobalMethods.texture_kg[row, col, layer, tex_class]; }
                                total_chem_weathered_mass_kg += weathered_mass_kg;
                                GlobalMethods.texture_kg[row, col, layer, tex_class] -= weathered_mass_kg;

                                //the following code accounts for the change in average size of the weathered class, 
                                //and the fact that a fraction of it therefore falls into a finer class as well

                                if (tex_class == 1)
                                {
                                    if (GlobalMethods.texture_kg[row, col, layer, tex_class] > 0.0000156252 * weathered_mass_kg)
                                    {

                                        GlobalMethods.texture_kg[row, col, layer, tex_class] -= 0.0000156252 * weathered_mass_kg;
                                        GlobalMethods.texture_kg[row, col, layer, tex_class + 1] += 0.0000156252 * weathered_mass_kg;
                                    }
                                    else
                                    {
                                        GlobalMethods.texture_kg[row, col, layer, tex_class + 1] += GlobalMethods.texture_kg[row, col, layer, tex_class];
                                        GlobalMethods.texture_kg[row, col, layer, tex_class] = 0;
                                    }
                                }
                                if (tex_class == 2)
                                {
                                    if (GlobalMethods.texture_kg[row, col, layer, tex_class] > 0.0000640041 * weathered_mass_kg)
                                    {
                                        GlobalMethods.texture_kg[row, col, layer, tex_class] -= 0.0000640041 * weathered_mass_kg;
                                        GlobalMethods.texture_kg[row, col, layer, tex_class + 1] += 0.0000640041 * weathered_mass_kg;
                                    }
                                    else
                                    {
                                        GlobalMethods.texture_kg[row, col, layer, tex_class + 1] += GlobalMethods.texture_kg[row, col, layer, tex_class];
                                        GlobalMethods.texture_kg[row, col, layer, tex_class] = 0;
                                    }
                                }
                                if (tex_class == 3)
                                {
                                    if (GlobalMethods.texture_kg[row, col, layer, tex_class] > 0.000125 * weathered_mass_kg)
                                    {
                                        GlobalMethods.texture_kg[row, col, layer, tex_class] -= 0.000125 * weathered_mass_kg;
                                        GlobalMethods.texture_kg[row, col, layer, tex_class + 1] += 0.000125 * weathered_mass_kg;
                                    }
                                    else
                                    {
                                        GlobalMethods.texture_kg[row, col, layer, tex_class + 1] += GlobalMethods.texture_kg[row, col, layer, tex_class];
                                        GlobalMethods.texture_kg[row, col, layer, tex_class] = 0;
                                    }
                                }
                                total_weath_mass += weathered_mass_kg;  //leached amount
                            }

                            //// clay neoformation
                            //fraction_neoform = neoform_constant * (Math.Exp(-Cfive * depth) - Math.Exp(-Csix * depth));
                            //if (fraction_neoform >= 1)
                            //{
                            //    Debug.WriteLine(" Warning: more than 100% of leached mass wants to become fine clay. This may indicate an error. Capping at 100%");
                            //    fraction_neoform = 1;
                            //}
                            //if (guiVariables.Daily_water) { fraction_neoform *= waterfactor[row,col]; }
                            //GlobalMethods.texture_kg[row, col, layer, 4] += total_weath_mass * fraction_neoform;
                            //total_fine_neoformed_mass_kg += total_weath_mass * fraction_neoform;
                            //total_weath_mass -= total_weath_mass * fraction_neoform;
                            //depth += GlobalMethods.layerthickness_m[row, col, layer] / 2;
                        }
                    }
                }
            }  //);
               //timeseries
            if (guiVariables.Timeseries.Total_chem_weath_checkbox)
            {
                timeseries_matrix[GlobalMethods.t, timeseries_order[16]] = total_chem_weathered_mass_kg;
            }
            if (guiVariables.Timeseries.Total_fine_formed_checkbox)
            {
                timeseries_matrix[GlobalMethods.t, timeseries_order[17]] = total_fine_neoformed_mass_kg;
            }

        }

        void soil_bioturbation()
        {

            Debug.WriteLine("\n--bioturbation--\n");
            try
            {
                //for bioturbation, we first calculate how much bioturbation (kg) this cell will experience, given its thickness
                //shallower soils do not experience the same amount as deeper soils
                //then we look at individual layers. Thicker layers, and layers closer to the surface, will experience more bioturbation kg
                //then, per layer, we will exchange the required bioturbation kg with the surrounding layers. 
                //Layers that are closer will exchange more than those that are further (regardless of whether they are deeper or closer to the surface)

                //XIA : For the Luxembourg study case, the kind of organic matter will affect the overall amount of bioturbation
                //Namely, 'young' organic matter (edible, from hornbeam) will allow more bioturbation than 'old' organic matter (inedible, from beech)
                //I simply calculate how much of each type the soil has, and then let the ratio between the two co-determine bioturbation.

                double local_bioturbation_kg, layer_bioturbation_kg, interlayer_bioturbation_kg;
                double layer_bio_activity_index, total_bio_activity_index, mass_distance_sum, mass_distance_layer;
                int layer, otherlayer;
                double fine_otherlayer_mass, fine_layer_mass;
                double total_soil_thickness_m;
                double depth, otherdepth, distance, potential_bioturbation_kg_m2_y;
                total_mass_bioturbed_kg = 0;
                double[,] temp_tex_som_kg = new double[GlobalMethods.max_soil_layers, 7]; // this will hold temporary changed values of texture until all bioturbation is done
                double[] layer_0 = new double[7], layer_0_after = new double[7];
                double mass_soil_before = 0, mass_soil_after = 0, mass_top_before = 0, mass_top_after = 0;
                // if (findnegativetexture()) { Debugger.Break(); }
                double lux_hornbeam_OM_litter_fraction = 0;

                double total_young_som_kg = 0, total_old_som_kg = 0;

                //another parallelization opportunity
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (GlobalMethods.t == 7 && row == 192 && col == 59) { diagnostic_mode = 1; }
                        else { diagnostic_mode = 0; }
                        if (GlobalMethods.dtm[row, col] != -9999 & GlobalMethods.soildepth_m[row, col] > 0)
                        {
                            remove_empty_layers(row, col);
                            update_all_soil_thicknesses(row, col);

                            mass_soil_before = total_soil_mass(row, col);
                            mass_top_before = total_layer_mass(row, col, 0);
                            total_soil_thickness_m = 0;
                            for (layer = 0; layer < GlobalMethods.max_soil_layers; layer++)
                            {
                                if (layer == 0)
                                {
                                    for (int tex = 0; tex < 5; tex++)
                                    {
                                        layer_0[tex] = GlobalMethods.texture_kg[row, col, 0, tex];
                                    }
                                    layer_0[5] = GlobalMethods.young_SOM_kg[row, col, 0];
                                    layer_0[6] = GlobalMethods.old_SOM_kg[row, col, 0];
                                }

                                //if (layer == 0 & !(GlobalMethods.layerthickness_m[row, col, layer] > 0)) { Debugger.Break(); }
                                //remove_empty_layers(row, col);
                                //if (layer == 0 & !(GlobalMethods.layerthickness_m[row, col, layer] > 0)) { Debugger.Break(); }
                                if (total_layer_mass(row, col, layer) > 0)  //this says: if the layer actually exists
                                {
                                    for (int prop = 0; prop < 5; prop++) { temp_tex_som_kg[layer, prop] = GlobalMethods.texture_kg[row, col, layer, prop]; }
                                    temp_tex_som_kg[layer, 5] = GlobalMethods.young_SOM_kg[row, col, layer];
                                    temp_tex_som_kg[layer, 6] = GlobalMethods.old_SOM_kg[row, col, layer];
                                    total_soil_thickness_m += GlobalMethods.layerthickness_m[row, col, layer];

                                    total_young_som_kg += GlobalMethods.young_SOM_kg[row, col, layer];
                                    total_old_som_kg += GlobalMethods.old_SOM_kg[row, col, layer];
                                }
                            } //can parallelize -PT

                            //LUX: in the lux case study, we need to know how much litter is hornbeam, i.e. palatable.
                            if (guiVariables.Version_lux_checkbox)
                            {
                                if (total_young_som_kg + total_old_som_kg > 0)
                                {
                                    lux_hornbeam_OM_litter_fraction = total_young_som_kg / (total_young_som_kg + total_old_som_kg);
                                }
                                else
                                {  // ArT quickfix attempt
                                    lux_hornbeam_OM_litter_fraction = 0.5;
                                }
                            }
                            //select vegetation parameters, same as GlobalMethods.creep
                            potential_bioturbation_kg_m2_y = 4.5;
                            if (guiVariables.Daily_water)
                            {
                                if (aridity_vegetation[row, col] < 1) { potential_bioturbation_kg_m2_y = 4 + 0.3; } // grassland
                                else { potential_bioturbation_kg_m2_y = 4 + 1.3; } // forest
                                                                                   // standard potential GlobalMethods.creep of 4 kg. 0.3 or 1.3 is added, based on vegetation type. Rates are derived from Wilkinson 2009: Breaking Ground and Gabet
                            }
                            // if (findnegativetexture()) { Debugger.Break(); }

                            // geen split voor voor depth decay voor verschillende vegetaties. Depth decay van GlobalMethods.creep aanhouden. 
                            //if(guiVariables.Daily_water)
                            //{
                            //    if (aridity_vegetation[row, col] < 1)
                            //    {
                            //        pot_bt_vegetation_kg = 0.3; // kg / m2/ y, from Gabet
                            //        depth_dec_vegetation = -1 / (0.5 / 2); // m-1, estimated root depth of 0.5 m
                            //    }
                            //    else
                            //    {
                            //        pot_bt_vegetation_kg = 1.3; // kg / m2/ y, from Gabet et al., 2003: https://doi.org/10.1146/annurev.earth.31.100901.141314
                            //        depth_dec_vegetation = -1 / (1.5 / 2); // m-1, estimated root depth of 1.5 m
                            //    }

                            //    pot_bt_animals_kg = 3; //  average animal burrowing rate, 30 ton/ha/yr, from Wilkinson et al., 2009: https://doi.org/10.1016/j.earscirev.2009.09.005 
                            //    depth_dec_animals = -1 / (1.0 / 2); // estimated, no source
                            //}
                            //// divide in animals and vegetation. animals is constant, vegetation differs per sort. Make new function?


                            //here we calculate the first quantity: how much bioturbation kg needs to happen in this location
                            local_bioturbation_kg = potential_bioturbation_kg_m2_y * (1 - Math.Exp(-bioturbation_depth_decay_constant * total_soil_thickness_m)) * GlobalMethods.dx * GlobalMethods.dx * dt;
                            if (local_bioturbation_kg < 0) // local_bt == 0 happens when soil is absent
                            {
                                Debug.WriteLine(" error in local_bioturbation calculation : zero mass");
                                Debug.WriteLine(" total soil thickness :" + total_soil_thickness_m + " at rc " + row + " " + col);
                                Debug.WriteLine("err_sbt1");

                            }
                            //LUX: if Luxembourg version: we assume that only hornbeam litter leads to bioturbation. More of it - more bioturbation.
                            if (guiVariables.Version_lux_checkbox) { local_bioturbation_kg *= lux_hornbeam_OM_litter_fraction; }

                            total_mass_bioturbed_kg += local_bioturbation_kg;



                            //now let's calculate layer-to-layer exchange to get to that local total needed.
                            depth = 0;
                            for (layer = 0; layer < GlobalMethods.max_soil_layers; layer++)
                            {

                                if (total_layer_fine_earth_mass(row, col, layer) > 0)  //this says: if there is actually fine earth in the layer. 
                                                                                       // That is necessary because we leave the stone fraction out of bioturbation
                                                                                       // and therefore purely stony layers are not involved in bioturbation
                                {
                                    //integration over the exponential decay function in JGR 2006 for the entire profile, and for the current layer.
                                    //then calculate the fraction of bioturbation that will happen in this layer, and multiply with total bioturbation in this cell
                                    fine_layer_mass = total_layer_fine_earth_mass(row, col, layer);
                                    layer_bio_activity_index = Math.Exp(-bioturbation_depth_decay_constant * depth) - (Math.Exp(-bioturbation_depth_decay_constant * (depth + GlobalMethods.layerthickness_m[row, col, layer])));
                                    total_bio_activity_index = 1 - (Math.Exp(-bioturbation_depth_decay_constant * total_soil_thickness_m));
                                    layer_bioturbation_kg = (layer_bio_activity_index / total_bio_activity_index) * local_bioturbation_kg;
                                    mass_distance_sum = 0;
                                    depth += GlobalMethods.layerthickness_m[row, col, layer] / 2;
                                    otherdepth = 0; distance = 0;

                                    if (GlobalMethods.layerthickness_m[row, col, layer] <= 0) { Debug.WriteLine(" error: layer thickness is 0 at GlobalMethods.t " + GlobalMethods.t + " r " + row + " c " + col); }

                                    //now that we know how much bioturbation originates in this layer,
                                    //let's look at other layers and decide which one of them exchanges how much of that good stuff.
                                    for (otherlayer = 0; otherlayer < GlobalMethods.max_soil_layers; otherlayer++)
                                    {
                                        if (total_layer_fine_earth_mass(row, col, otherlayer) > 0)  //this says: if there is actually fine earth in the layer.
                                                                                                    // That is necessary because we leave the stone fraction out of bioturbation
                                                                                                    // and therefore purely stony layers are not involved in bioturbation
                                        {

                                            otherdepth += GlobalMethods.layerthickness_m[row, col, otherlayer] / 2;
                                            distance = Math.Abs(otherdepth - depth);

                                            if (distance < 0) { Debug.WriteLine(" distance between layers is 0 m at row " + row + " col " + col + " layerdepth " + depth + " otherlayerdepth " + otherdepth); }

                                            if (double.IsNaN(GlobalMethods.texture_kg[row, col, otherlayer, 1])) { Debug.WriteLine(" texture 1 NaN " + GlobalMethods.t + " rc " + row + "  " + col + " layers " + layer + " " + otherlayer); }
                                            if (double.IsNaN(GlobalMethods.texture_kg[row, col, otherlayer, 2])) { Debug.WriteLine(" texture 2 NaN " + GlobalMethods.t + " rc " + row + "  " + col + " layers " + layer + " " + otherlayer); }
                                            if (double.IsNaN(GlobalMethods.texture_kg[row, col, otherlayer, 3])) { Debug.WriteLine(" texture 3 NaN " + GlobalMethods.t + " rc " + row + "  " + col + " layers " + layer + " " + otherlayer); }
                                            if (double.IsNaN(GlobalMethods.texture_kg[row, col, otherlayer, 4])) { Debug.WriteLine(" texture 4 NaN " + GlobalMethods.t + " rc " + row + "  " + col + " layers " + layer + " " + otherlayer); }
                                            if (double.IsNaN(GlobalMethods.young_SOM_kg[row, col, otherlayer])) { Debug.WriteLine(" young som NaN " + GlobalMethods.t + " rc " + row + "  " + col + " layers " + layer + " " + otherlayer); }
                                            if (double.IsNaN(GlobalMethods.old_SOM_kg[row, col, otherlayer])) { Debug.WriteLine(" old som NaN " + GlobalMethods.t + " rc " + row + "  " + col + " layers " + layer + " " + otherlayer); }

                                            if ((GlobalMethods.texture_kg[row, col, otherlayer, 1] < 0)) { Debug.WriteLine(" texture 1 null " + GlobalMethods.t + " rc " + row + "  " + col + " layers " + layer + " " + otherlayer); }
                                            if ((GlobalMethods.texture_kg[row, col, otherlayer, 2] < 0)) { Debug.WriteLine(" texture 2 null " + GlobalMethods.t + " rc " + row + "  " + col + " layers " + layer + " " + otherlayer); }
                                            if ((GlobalMethods.texture_kg[row, col, otherlayer, 3] < 0)) { Debug.WriteLine(" texture 3 null " + GlobalMethods.t + " rc " + row + "  " + col + " layers " + layer + " " + otherlayer); }
                                            if ((GlobalMethods.texture_kg[row, col, otherlayer, 4] < 0)) { Debug.WriteLine(" texture 4 null " + GlobalMethods.t + " rc " + row + "  " + col + " layers " + layer + " " + otherlayer); }
                                            if ((GlobalMethods.young_SOM_kg[row, col, otherlayer] < 0)) { Debug.WriteLine(" young som null " + GlobalMethods.t + " rc " + row + "  " + col + " layers " + layer + " " + otherlayer); }
                                            if ((GlobalMethods.old_SOM_kg[row, col, otherlayer] < 0)) { Debug.WriteLine(" old som null " + GlobalMethods.t + " rc " + row + "  " + col + " layers " + layer + " " + otherlayer); }

                                            if (otherlayer != layer) { mass_distance_sum += (GlobalMethods.texture_kg[row, col, otherlayer, 1] + GlobalMethods.texture_kg[row, col, otherlayer, 2] + GlobalMethods.texture_kg[row, col, otherlayer, 3] + GlobalMethods.texture_kg[row, col, otherlayer, 4] + GlobalMethods.young_SOM_kg[row, col, otherlayer] + GlobalMethods.old_SOM_kg[row, col, otherlayer]) / distance; }

                                            otherdepth += GlobalMethods.layerthickness_m[row, col, otherlayer] / 2;
                                            if (double.IsNaN(mass_distance_sum)) { Debug.WriteLine(" B NaN mass distance in bioturbation GlobalMethods.t " + GlobalMethods.t + " rc " + row + "  " + col + " layers " + layer + " " + otherlayer); }
                                            if (double.IsNaN(distance)) { Debug.WriteLine(" NaN  distance in bioturbation GlobalMethods.t " + GlobalMethods.t + " rc " + row + "  " + col + " layers " + layer + " " + otherlayer); }
                                            if (double.IsNaN((GlobalMethods.layerthickness_m[row, col, otherlayer] / 2))) { Debug.WriteLine(" NaN layerthick in bioturbation GlobalMethods.t " + GlobalMethods.t + " rc " + row + "  " + col + " layers " + layer + " " + otherlayer); }
                                        }
                                    }
                                    double check_mass_distance = 0;
                                    otherdepth = 0; distance = 0;
                                    double BT_fraction = 0;
                                    for (otherlayer = 0; otherlayer < GlobalMethods.max_soil_layers; otherlayer++)
                                    {
                                        otherdepth += GlobalMethods.layerthickness_m[row, col, otherlayer] / 2;
                                        if (total_layer_fine_earth_mass(row, col, otherlayer) > 0 && layer != otherlayer)  //this says: if the other layer actually exists and if it's not the current layer
                                        {

                                            distance = Math.Abs(otherdepth - depth);
                                            if (layer != otherlayer)
                                            {
                                                mass_distance_layer = (GlobalMethods.texture_kg[row, col, otherlayer, 1] + GlobalMethods.texture_kg[row, col, otherlayer, 2] + GlobalMethods.texture_kg[row, col, otherlayer, 3] + GlobalMethods.texture_kg[row, col, otherlayer, 4] + GlobalMethods.young_SOM_kg[row, col, otherlayer] + GlobalMethods.old_SOM_kg[row, col, otherlayer]) / distance;
                                            }
                                            else
                                            {
                                                mass_distance_layer = (GlobalMethods.texture_kg[row, col, otherlayer, 1] + GlobalMethods.texture_kg[row, col, otherlayer, 2] + GlobalMethods.texture_kg[row, col, otherlayer, 3] + GlobalMethods.texture_kg[row, col, otherlayer, 4] + GlobalMethods.young_SOM_kg[row, col, otherlayer] + GlobalMethods.old_SOM_kg[row, col, otherlayer]) / (GlobalMethods.layerthickness_m[row, col, otherlayer] / 2);
                                            }
                                            if (double.IsNaN(mass_distance_layer))
                                            {
                                                Debug.WriteLine(" NaN mass distance layer in bioturbation GlobalMethods.t " + GlobalMethods.t + " rc " + row + "  " + col + " layers " + layer + " " + otherlayer); Debug.WriteLine("err_sbt2");
                                            }
                                            if (mass_distance_sum == 0) { Debug.WriteLine(" zero mass distance sum"); }
                                            //here we calculate the amount of material bioturbated between the current layer and the current otherlayer
                                            interlayer_bioturbation_kg = layer_bioturbation_kg * (mass_distance_layer / mass_distance_sum);
                                            check_mass_distance += mass_distance_layer / mass_distance_sum;
                                            BT_fraction += mass_distance_layer / mass_distance_sum;
                                            if (interlayer_bioturbation_kg < 0)
                                            {
                                                Debug.WriteLine("err_sbt3");
                                            }
                                            if (double.IsNaN(interlayer_bioturbation_kg))
                                            {
                                                Debug.WriteLine(" NaN interlayer bioturbation kg in bioturbation GlobalMethods.t " + GlobalMethods.t + " rc " + row + " " + col + " layers " + layer + " " + otherlayer);
                                                Debug.WriteLine(" " + mass_distance_layer + " " + mass_distance_sum + " " + layer_bioturbation_kg);
                                                Debug.WriteLine("err_sbt4");

                                            }
                                            fine_otherlayer_mass = GlobalMethods.texture_kg[row, col, otherlayer, 1] + GlobalMethods.texture_kg[row, col, otherlayer, 2] + GlobalMethods.texture_kg[row, col, otherlayer, 3] + GlobalMethods.texture_kg[row, col, otherlayer, 4] + GlobalMethods.young_SOM_kg[row, col, otherlayer] + GlobalMethods.old_SOM_kg[row, col, otherlayer];
                                            if (double.IsNaN(fine_otherlayer_mass)) { Debug.WriteLine(" NaN fine otherlayer mass in bioturbation "); }
                                            if ((fine_otherlayer_mass <= 0))
                                            {
                                                Debug.WriteLine(" fineotherlayermass " + fine_otherlayer_mass + " GlobalMethods.t " + GlobalMethods.t + " rc " + row + "  " + col + " layers " + layer + " " + otherlayer);
                                                Debug.WriteLine("err_sbt5");
                                            }

                                            //weathered_mass_kg may be more than present in the other layer, the current layer, or both - in that case one or both of the layers will become mixtures of the original two layers
                                            double fromlayertomixture_kg = 0, fromotherlayertomixture_kg = 0, totalmixturemass_kg = 0, massfromlayer = 0, massfromotherlayer = 0, dmass_l, dmass_ol;
                                            double[] mixture_kg = new double[7];
                                            fromlayertomixture_kg = Math.Min(fine_layer_mass, (interlayer_bioturbation_kg / 2));
                                            fromotherlayertomixture_kg = Math.Min(fine_otherlayer_mass, (interlayer_bioturbation_kg / 2));
                                            // totalmixturemass_kg = fromlayertomixture_kg + fromotherlayertomixture_kg;

                                            if ((fromlayertomixture_kg + fromotherlayertomixture_kg) > 1E-6)  // if there is actual exchange (which is not the case when all fine material is removed)
                                            {
                                                //now add to mixture, and take from donors
                                                double massin_l = 0, massin_ol = 0;
                                                // texture
                                                for (int prop = 1; prop < 5; prop++)
                                                {
                                                    //checks
                                                    if (temp_tex_som_kg[layer, prop] < 0)
                                                    {
                                                        Debug.WriteLine("err_sbt6");
                                                    }
                                                    if (temp_tex_som_kg[otherlayer, prop] < 0)
                                                    {
                                                        Debug.WriteLine("err_sbt7");
                                                    }

                                                    //determine how much mass can be exchanged,. Do not take more than is present in the temporary layer to prevent negative textures in the end
                                                    //Should not happen, mass of top layer should stay constant, but happens anyway
                                                    dmass_l = (fromlayertomixture_kg / fine_layer_mass) * GlobalMethods.texture_kg[row, col, layer, prop];
                                                    dmass_ol = (fromotherlayertomixture_kg / fine_otherlayer_mass) * GlobalMethods.texture_kg[row, col, otherlayer, prop];

                                                    if (dmass_l > temp_tex_som_kg[layer, prop]) { dmass_l = temp_tex_som_kg[layer, prop]; }
                                                    if (dmass_ol > temp_tex_som_kg[otherlayer, prop]) { dmass_ol = temp_tex_som_kg[otherlayer, prop]; }

                                                    //take mass from donors to mix
                                                    mixture_kg[prop] += (dmass_l + dmass_ol);
                                                    massfromlayer += dmass_l;
                                                    massfromotherlayer += dmass_ol;

                                                    temp_tex_som_kg[layer, prop] -= dmass_l;
                                                    temp_tex_som_kg[otherlayer, prop] -= dmass_ol;

                                                    //mixture_kg[prop] += (fromlayertomixture_kg / fine_layer_mass) * GlobalMethods.texture_kg[row, col, layer, prop];
                                                    //if ((fromlayertomixture_kg / fine_layer_mass) * GlobalMethods.texture_kg[row, col, layer, prop] < 0) { Debugger.Break(); }
                                                    //massfromlayer += (fromlayertomixture_kg / fine_layer_mass) * GlobalMethods.texture_kg[row, col, layer, prop];
                                                    //mixture_kg[prop] += (fromotherlayertomixture_kg / fine_otherlayer_mass) * GlobalMethods.texture_kg[row, col, otherlayer, prop];
                                                    //massfromotherlayer += (fromotherlayertomixture_kg / fine_otherlayer_mass) * GlobalMethods.texture_kg[row, col, otherlayer, prop];
                                                    //if ((fromotherlayertomixture_kg / fine_otherlayer_mass) * GlobalMethods.texture_kg[row, col, otherlayer, prop] < 0) { Debugger.Break(); }

                                                    //temp_tex_som_kg[otherlayer, prop] -= (fromotherlayertomixture_kg / fine_otherlayer_mass) * GlobalMethods.texture_kg[row, col, otherlayer, prop];
                                                    //temp_tex_som_kg[layer, prop] -= (fromlayertomixture_kg / fine_layer_mass) * GlobalMethods.texture_kg[row, col, layer, prop];

                                                    if (temp_tex_som_kg[layer, prop] < 0)
                                                    {
                                                        Debug.WriteLine("err_sbt8");
                                                    }
                                                    if (temp_tex_som_kg[otherlayer, prop] < 0)
                                                    {
                                                        Debug.WriteLine("err_sbt9");
                                                    }

                                                }
                                                //young OM
                                                dmass_l = (fromlayertomixture_kg / fine_layer_mass) * (GlobalMethods.young_SOM_kg[row, col, layer]);
                                                dmass_ol = (fromotherlayertomixture_kg / fine_otherlayer_mass) * (GlobalMethods.young_SOM_kg[row, col, otherlayer]);

                                                if (dmass_l > temp_tex_som_kg[layer, 5]) { dmass_l = temp_tex_som_kg[layer, 5]; }
                                                if (dmass_ol > temp_tex_som_kg[otherlayer, 5]) { dmass_ol = temp_tex_som_kg[otherlayer, 5]; }

                                                //take mass from donors to mix
                                                mixture_kg[5] += (dmass_l + dmass_ol);
                                                massfromlayer += dmass_l;
                                                massfromotherlayer += dmass_ol;

                                                temp_tex_som_kg[layer, 5] -= dmass_l;
                                                temp_tex_som_kg[otherlayer, 5] -= dmass_ol;

                                                //old OM
                                                // if (layer == 0) { Debugger.Break(); }
                                                dmass_l = (fromlayertomixture_kg / fine_layer_mass) * (GlobalMethods.old_SOM_kg[row, col, layer]);
                                                dmass_ol = (fromotherlayertomixture_kg / fine_otherlayer_mass) * (GlobalMethods.old_SOM_kg[row, col, otherlayer]);

                                                if (dmass_l > temp_tex_som_kg[layer, 6]) { dmass_l = temp_tex_som_kg[layer, 6]; }
                                                if (dmass_ol > temp_tex_som_kg[otherlayer, 6]) { dmass_ol = temp_tex_som_kg[otherlayer, 6]; }

                                                //take mass from donors to mix
                                                mixture_kg[6] += (dmass_l + dmass_ol);
                                                massfromlayer += dmass_l;
                                                massfromotherlayer += dmass_ol;

                                                temp_tex_som_kg[layer, 6] -= dmass_l;
                                                temp_tex_som_kg[otherlayer, 6] -= dmass_ol;

                                                // checks
                                                if (temp_tex_som_kg[layer, 5] < 0)
                                                {
                                                    Debug.WriteLine("err_sbt10");
                                                }
                                                if (temp_tex_som_kg[otherlayer, 5] < 0)
                                                {
                                                    Debug.WriteLine("err_sbt11");
                                                }

                                                //now give from mixture to receivers
                                                totalmixturemass_kg = massfromlayer + massfromotherlayer;
                                                if (totalmixturemass_kg == 0)
                                                {
                                                    Debug.WriteLine("err_sbt11");
                                                }

                                                // if (findnegativetexture()) { Debugger.Break(); }


                                                for (int prop = 1; prop < 7; prop++)
                                                {
                                                    if (temp_tex_som_kg[layer, prop] < 0)
                                                    {
                                                        Debug.WriteLine("err_sbt12");
                                                    }
                                                    if (temp_tex_som_kg[otherlayer, prop] < 0)
                                                    {
                                                        Debug.WriteLine("err_sbt13");
                                                    }

                                                    if (mixture_kg[prop] < 0)
                                                    {
                                                        Debug.WriteLine("err_sbt14");
                                                    }
                                                    temp_tex_som_kg[otherlayer, prop] += mixture_kg[prop] * (massfromotherlayer / totalmixturemass_kg);
                                                    massin_ol += mixture_kg[prop] * (massfromotherlayer / totalmixturemass_kg);
                                                    temp_tex_som_kg[layer, prop] += mixture_kg[prop] * (massfromlayer / totalmixturemass_kg);
                                                    massin_l += mixture_kg[prop] * (massfromlayer / totalmixturemass_kg);
                                                    //mixture_kg[prop] = 0;  // that's not perse needed, but feels clean

                                                    if (temp_tex_som_kg[layer, prop] < 0)
                                                    {
                                                        Debug.WriteLine("err_sbt15");
                                                    }
                                                    if (temp_tex_som_kg[otherlayer, prop] < 0)
                                                    {
                                                        Debug.WriteLine("err_sbt16");
                                                    }
                                                }
                                            }


                                            //all sorts of checks - we should never have values under zero, or NotANumber NaN
                                            if (temp_tex_som_kg[otherlayer, 1] < 0)
                                            {
                                                Debug.WriteLine(" texture 1 null " + GlobalMethods.t + " rc " + row + "  " + col + " otherlayers " + layer + " (" + total_layer_mass(row, col, layer) + "kg) " + otherlayer + " (" + total_layer_mass(row, col, otherlayer) + "kg) ");
                                            }
                                            if (temp_tex_som_kg[otherlayer, 2] < 0) { Debug.WriteLine(" texture 2 null " + GlobalMethods.t + " rc " + row + "  " + col + " otherlayers " + layer + " " + otherlayer); }
                                            if (temp_tex_som_kg[otherlayer, 3] < 0) { Debug.WriteLine(" texture 3 null " + GlobalMethods.t + " rc " + row + "  " + col + " otherlayers " + layer + " " + otherlayer); }
                                            if (temp_tex_som_kg[otherlayer, 4] < 0) { Debug.WriteLine(" texture 4 null " + GlobalMethods.t + " rc " + row + "  " + col + " otherlayers " + layer + " " + otherlayer); }
                                            if (temp_tex_som_kg[otherlayer, 5] < 0) { Debug.WriteLine(" young som null " + GlobalMethods.t + " rc " + row + "  " + col + " otherlayers " + layer + " " + otherlayer); }
                                            if (temp_tex_som_kg[otherlayer, 6] < 0) { Debug.WriteLine(" old som null " + GlobalMethods.t + " rc " + row + "  " + col + " otherlayers " + layer + " " + otherlayer); }

                                            if (temp_tex_som_kg[layer, 1] < 0) { Debug.WriteLine(" texture 1 null " + GlobalMethods.t + " rc " + row + "  " + col + " layer " + layer + " " + otherlayer); }
                                            if (temp_tex_som_kg[layer, 2] < 0) { Debug.WriteLine(" texture 2 null " + GlobalMethods.t + " rc " + row + "  " + col + " layer " + layer + " " + otherlayer); }
                                            if (temp_tex_som_kg[layer, 3] < 0) { Debug.WriteLine(" texture 3 null " + GlobalMethods.t + " rc " + row + "  " + col + " layer " + layer + " " + otherlayer); }
                                            if (temp_tex_som_kg[layer, 4] < 0) { Debug.WriteLine(" texture 4 null " + GlobalMethods.t + " rc " + row + "  " + col + " layer " + layer + " " + otherlayer); }
                                            if (temp_tex_som_kg[layer, 5] < 0) { Debug.WriteLine(" young som null " + GlobalMethods.t + " rc " + row + "  " + col + " layer " + layer + " " + otherlayer); }
                                            if (temp_tex_som_kg[layer, 6] < 0) { Debug.WriteLine(" old som null " + GlobalMethods.t + " rc " + row + "  " + col + " layer " + layer + " " + otherlayer); }

                                            if (double.IsNaN(temp_tex_som_kg[otherlayer, 1]))
                                            {
                                                Debug.WriteLine(" texture 1 NaN " + GlobalMethods.t + " rc " + row + "  " + col + " otherlayers " + layer + " (" + total_layer_mass(row, col, layer) + "kg) " + otherlayer + " (" + total_layer_mass(row, col, otherlayer) + "kg) ");
                                            }
                                            if (double.IsNaN(temp_tex_som_kg[otherlayer, 2])) { Debug.WriteLine(" texture 2 NaN " + GlobalMethods.t + " rc " + row + "  " + col + " otherlayers " + layer + " " + otherlayer); }
                                            if (double.IsNaN(temp_tex_som_kg[otherlayer, 3])) { Debug.WriteLine(" texture 3 NaN " + GlobalMethods.t + " rc " + row + "  " + col + " otherlayers " + layer + " " + otherlayer); }
                                            if (double.IsNaN(temp_tex_som_kg[otherlayer, 4])) { Debug.WriteLine(" texture 4 NaN " + GlobalMethods.t + " rc " + row + "  " + col + " otherlayers " + layer + " " + otherlayer); }
                                            if (double.IsNaN(temp_tex_som_kg[otherlayer, 5])) { Debug.WriteLine(" young som NaN " + GlobalMethods.t + " rc " + row + "  " + col + " otherlayers " + layer + " " + otherlayer); }
                                            if (double.IsNaN(temp_tex_som_kg[otherlayer, 6])) { Debug.WriteLine(" old som NaN " + GlobalMethods.t + " rc " + row + "  " + col + " otherlayers " + layer + " " + otherlayer); }

                                            if (double.IsNaN(temp_tex_som_kg[layer, 1])) { Debug.WriteLine(" texture 1 NaN " + GlobalMethods.t + " rc " + row + "  " + col + " layer " + layer + " " + otherlayer); }
                                            if (double.IsNaN(temp_tex_som_kg[layer, 2])) { Debug.WriteLine(" texture 2 NaN " + GlobalMethods.t + " rc " + row + "  " + col + " layer " + layer + " " + otherlayer); }
                                            if (double.IsNaN(temp_tex_som_kg[layer, 3])) { Debug.WriteLine(" texture 3 NaN " + GlobalMethods.t + " rc " + row + "  " + col + " layer " + layer + " " + otherlayer); }
                                            if (double.IsNaN(temp_tex_som_kg[layer, 4])) { Debug.WriteLine(" texture 4 NaN " + GlobalMethods.t + " rc " + row + "  " + col + " layer " + layer + " " + otherlayer); }
                                            if (double.IsNaN(temp_tex_som_kg[layer, 5])) { Debug.WriteLine(" young som NaN " + GlobalMethods.t + " rc " + row + "  " + col + " layer " + layer + " " + otherlayer); }
                                            if (double.IsNaN(temp_tex_som_kg[layer, 6])) { Debug.WriteLine(" old som NaN " + GlobalMethods.t + " rc " + row + "  " + col + " layer " + layer + " " + otherlayer); }
                                        }
                                        otherdepth += GlobalMethods.layerthickness_m[row, col, otherlayer] / 2;
                                    }
                                    //if (Math.Round(check_mass_distance,4) != 1) { Debugger.Break(); }
                                    //if (findnegativetexture()) { Debugger.Break(); }

                                    // if (Math.Round(BT_fraction, 6) != 1) { Debugger.Break(); }
                                    depth += GlobalMethods.layerthickness_m[row, col, layer] / 2;
                                }
                            } // end for layer
                              // now we know the new, bioturbated amounts in every layer in this row col, let's store them in the main GlobalMethods.texture_kg variables
                            for (layer = 0; layer < GlobalMethods.max_soil_layers; layer++)
                            {
                                if (layer == 0 & temp_tex_som_kg[0, 2] == 0)
                                {
                                    //Debug.WriteLine("err_sbt_16a. Possible empty top layer after BT 0: {0}, {1}, {2}, {3}, {4}, {5}, {6}. GlobalMethods.t {7}, row {8}, col {9}, dlayer {10}", layer_0[0], layer_0[1], layer_0[2], layer_0[3], layer_0[4], layer_0[5], layer_0[6], GlobalMethods.t, row, col, GlobalMethods.layerthickness_m[row, col, 0]);
                                    //this does not really test for an empty top layer - just for a silt-less top layer.
                                    if (layer == 0 & total_layer_fine_earth_mass(row, col, 0) == 0)
                                    {
                                        //this does!
                                        //Debug.WriteLine("confirmed!");
                                    }
                                }
                                for (int prop = 1; prop < 5; prop++)
                                {
                                    if (temp_tex_som_kg[layer, prop] < 0)
                                    {
                                        Debug.WriteLine("err_sbt17");
                                    }
                                    GlobalMethods.texture_kg[row, col, layer, prop] = temp_tex_som_kg[layer, prop];
                                    layer_0_after[prop] = temp_tex_som_kg[layer, prop];
                                    temp_tex_som_kg[layer, prop] = 0;
                                }
                                GlobalMethods.young_SOM_kg[row, col, layer] = temp_tex_som_kg[layer, 5];
                                GlobalMethods.old_SOM_kg[row, col, layer] = temp_tex_som_kg[layer, 6];
                                layer_0_after[5] = temp_tex_som_kg[layer, 5];
                                layer_0_after[6] = temp_tex_som_kg[layer, 6];
                                temp_tex_som_kg[layer, 5] = 0;
                                temp_tex_som_kg[layer, 6] = 0;
                            } //end for layer
                              // if (findnegativetexture()) { Debugger.Break(); }

                            mass_soil_after = total_soil_mass(row, col);
                            mass_top_after = total_layer_mass(row, col, 0);

                            if (Math.Abs(mass_soil_before - mass_soil_after) > 1E-8 | Math.Abs(mass_top_before - mass_top_after) > 1E-8)
                            {
                                Debug.WriteLine("Mass loss during bioturbation");
                                // Debugger.Break(); 
                            }

                        } // end GlobalMethods.dtm!=-9999
                    }// for col
                } // end for row
                  // if (findnegativetexture()) { Debugger.Break(); }


                if (guiVariables.Timeseries.Total_mass_bioturbed_checkbox)
                {
                    timeseries_matrix[GlobalMethods.t, timeseries_order[19]] = total_mass_bioturbed_kg;
                }
                if (NA_in_map(GlobalMethods.dtm) > 0 | NA_in_map(GlobalMethods.soildepth_m) > 0)
                {
                    Debug.WriteLine("err_sbt20");
                }

            }
            catch { Debug.WriteLine(" No valid text in textbox bioturbation "); }

        } // nieuwe code van Arnaud

        void soil_litter_cycle()
        {
            // uses parameters from Carbon Cycle for now
            try
            {
                double litter_input_kg;

                //this line keeps young (hornbeam) OM completely gone from the surface every second year (reflecting that,
                //in reality, part of the year is unprotected). MvdM I added the else to rest the decomposition rate
                if (GlobalMethods.t % 2 == 0) { potential_young_decomp_rate = 1; } else { potential_young_decomp_rate = Convert.ToDouble(carbon_y_decomp_rate_textbox.Text); }

                calculate_TPI(7);
                double a = -0.33;
                double b = 28.33;
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        //calculating hornbeam fraction
                        GlobalMethods.hornbeam_cover_fraction[row, col] = 1 - Math.Exp(a + b * GlobalMethods.tpi[row, col]) / (1 + Math.Exp(a + b * GlobalMethods.tpi[row, col]));


                        litter_input_kg = potential_OM_input; // MvdM no changes in litter input due to soil thickness. 

                        // All litter to litter layer matrix
                        GlobalMethods.litter_kg[row, col, 0] += .230 * (GlobalMethods.hornbeam_cover_fraction[row, col]); // Hornbeam
                        GlobalMethods.litter_kg[row, col, 1] += .403 * (1 - GlobalMethods.hornbeam_cover_fraction[row, col]); // Beech

                        // Decomposition
                        GlobalMethods.litter_kg[row, col, 0] *= (1 - .97); // Hornbeam
                        GlobalMethods.litter_kg[row, col, 1] *= (1 - .45); // Beech
                    }
                }
            }
            catch { Debug.WriteLine(" Crash in litter cycle "); }
        }


        void soil_carbon_cycle()
        {
            try
            {
                double local_OM_input_kg, layer_OM_input_kg;
                double young_decomposition_rate, old_decomposition_rate;
                //Debug.WriteLine("succesfully read parameters for soil carbon");
                double depth;
                double total_soil_thickness;
                double layer_OM_input_index, total_OM_input_index;
                double dz_dx, dz_dy, slope_local, dynamic_TWI;
                int layer;
                total_OM_input_kg = 0;
                if (guiVariables.Version_lux_checkbox)
                {
                    calculate_TPI(7);
                    double a = -0.33;
                    double b = 28.33;
                    for (int row = 0; row < GlobalMethods.nr; row++)
                    {
                        for (int col = 0; col < GlobalMethods.nc; col++)
                        {
                            //calculating hornbeam fraction
                            GlobalMethods.hornbeam_cover_fraction[row, col] = 1 - Math.Exp(a + b * GlobalMethods.tpi[row, col]) / (1 + Math.Exp(a + b * GlobalMethods.tpi[row, col]));
                        }
                    }
                    //this line keeps young (hornbeam) OM completely gone from the surface every second year (reflecting that,
                    //in reality, part of the year is unprotected). MvdM I added the else to rest the decomposition rate
                    if (GlobalMethods.t % 2 == 0) { potential_young_decomp_rate = 1; } else { potential_young_decomp_rate = Convert.ToDouble(carbon_y_decomp_rate_textbox.Text); }
                }

                if (NA_in_map(GlobalMethods.dtm) > 0 | NA_in_map(GlobalMethods.soildepth_m) > 0)
                {
                    Debug.WriteLine("err_cc1");
                }
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    //Parallel.For(0, GlobalMethods.nc, i =>                    //we parallelize over cols
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (guiVariables.Daily_water)
                        {
                            if (aridity_vegetation[row, col] < 1) { potential_OM_input = 0.67; } // grassland
                            else { potential_OM_input = 0.62; } // forest
                            if (GlobalMethods.t > (guiVariables.End_time - 500)) { potential_OM_input = 0.42; } // arable land
                        }
                        total_soil_thickness = 0;
                        for (layer = 0; layer < GlobalMethods.max_soil_layers; layer++)
                        {
                            if (GlobalMethods.layerthickness_m[row, col, layer] > 0)
                            {
                                total_soil_thickness += GlobalMethods.layerthickness_m[row, col, layer];
                            }
                        }
                        local_OM_input_kg = potential_OM_input * (1 - Math.Exp(-OM_input_depth_decay_constant * total_soil_thickness)) * GlobalMethods.dx * GlobalMethods.dx * dt;
                        total_OM_input_kg += local_OM_input_kg;
                        depth = 0;


                        for (layer = 0; layer < GlobalMethods.max_soil_layers; layer++)
                        {
                            if (GlobalMethods.layerthickness_m[row, col, layer] > 0)
                            {
                                // if (layer == 0) { Debugger.Break(); }
                                layer_OM_input_index = -1 / OM_input_depth_decay_constant * (Math.Exp(-OM_input_depth_decay_constant * (depth + GlobalMethods.layerthickness_m[row, col, layer])) - Math.Exp(-OM_input_depth_decay_constant * depth));
                                total_OM_input_index = -1 / OM_input_depth_decay_constant * (Math.Exp(-OM_input_depth_decay_constant * (total_soil_thickness)) - 1);
                                layer_OM_input_kg = (layer_OM_input_index / total_OM_input_index) * local_OM_input_kg;

                                if (guiVariables.Version_lux_checkbox)
                                {
                                    GlobalMethods.young_SOM_kg[row, col, layer] += layer_OM_input_kg * (GlobalMethods.hornbeam_cover_fraction[row, col]);
                                    GlobalMethods.old_SOM_kg[row, col, layer] += layer_OM_input_kg * (1 - GlobalMethods.hornbeam_cover_fraction[row, col]);
                                }
                                else
                                {
                                    GlobalMethods.young_SOM_kg[row, col, layer] += layer_OM_input_kg * (1 - humification_fraction);
                                    GlobalMethods.old_SOM_kg[row, col, layer] += layer_OM_input_kg * (humification_fraction);
                                }
                                if (double.IsNaN(GlobalMethods.young_SOM_kg[row, col, layer]))
                                {
                                    Debug.WriteLine("err_cc2");
                                }
                                //decomposition gets lost as CO2 to the air (and soil water)
                                depth += GlobalMethods.layerthickness_m[row, col, layer] / 2;
                                //young_decomposition_rate = potential_young_decomp_rate * Math.Exp(-young_CTI_decay_constant * dynamic_TWI) * Math.Exp(-young_depth_decay_constant * depth);
                                //old_decomposition_rate = potential_old_decomp_rate * Math.Exp(-old_CTI_decay_constant * dynamic_TWI) * Math.Exp(-old_depth_decay_constant * depth);
                                young_decomposition_rate = potential_young_decomp_rate * 1 * Math.Exp(-young_depth_decay_constant * depth);
                                old_decomposition_rate = potential_old_decomp_rate * 1 * Math.Exp(-old_depth_decay_constant * depth);
                                GlobalMethods.young_SOM_kg[row, col, layer] *= (1 - young_decomposition_rate);
                                GlobalMethods.old_SOM_kg[row, col, layer] *= (1 - old_decomposition_rate);
                                //Debug.WriteLine(" cell  " + row + " " + col + " has layer_OM_input of " + layer_OM_input_kg);
                                depth += GlobalMethods.layerthickness_m[row, col, layer] / 2;
                                if (GlobalMethods.young_SOM_kg[row, col, layer] < 0 | GlobalMethods.old_SOM_kg[row, col, layer] < 0)
                                {
                                    Debug.WriteLine("err_cc3");
                                }
                            }

                        }
                    }

                }
                if (guiVariables.Timeseries.Total_OM_input_checkbox)
                {
                    timeseries_matrix[GlobalMethods.t, timeseries_order[20]] = total_OM_input_kg;
                }
                if (NA_in_map(GlobalMethods.dtm) > 0 | NA_in_map(GlobalMethods.soildepth_m) > 0)
                {
                    Debug.WriteLine("err_cc4");
                }

            }
            catch { Debug.WriteLine(" Crash in soil carbon cycle "); }

        }

        void soil_clay_translocation()
        {
            //possibly a function of local wetness / infiltration, but for now not/.
            double Iavg = 0, Imin = 10000000, Imax = 0;

            if (guiVariables.Daily_water)
            {
                int Icount = 0;
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (GlobalMethods.dtm[row, col] != -9999)
                        {
                            if (Imin > Iy[row, col]) { Imin = Iy[row, col]; }
                            if (Imax < Iy[row, col]) { Imax = Iy[row, col]; }
                            Iavg += Iy[row, col];
                            Icount++;
                        }
                    }
                }
                Iavg /= Icount;
            }

            int layer;
            double eluviated_kg, depth;
            total_fine_eluviated_mass_kg = 0;
            try
            {
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {


                        depth = 0;
                        for (layer = 0; layer < GlobalMethods.max_soil_layers - 1; layer++)   // we loop through all layers except the lower one - clay translocation there has no lower recipient
                        {
                            if (GlobalMethods.layerthickness_m[row, col, layer] > 0 && GlobalMethods.layerthickness_m[row, col, layer + 1] > 0)  // both source and sink layers have to exist.
                            {
                                if (GlobalMethods.texture_kg[row, col, layer, 4] > 0)
                                {
                                    depth += GlobalMethods.layerthickness_m[row, col, layer] / 2;
                                    double totalweight = GlobalMethods.texture_kg[row, col, layer, 0] + GlobalMethods.texture_kg[row, col, layer, 1] + GlobalMethods.texture_kg[row, col, layer, 2] + GlobalMethods.texture_kg[row, col, layer, 3] + GlobalMethods.texture_kg[row, col, layer, 4] + GlobalMethods.young_SOM_kg[row, col, layer] + GlobalMethods.old_SOM_kg[row, col, layer];
                                    //calculate the mass of eluviation
                                    if (CT_depth_decay_checkbox.Checked) { eluviated_kg = max_eluviation * (1 - Math.Exp(-Cclay * GlobalMethods.texture_kg[row, col, layer, 4] / totalweight)) * Math.Exp(-ct_depthdec * depth) * dt * GlobalMethods.dx * GlobalMethods.dx; }
                                    else { eluviated_kg = max_eluviation * (1 - Math.Exp(-Cclay * GlobalMethods.texture_kg[row, col, layer, 4] / totalweight)) * dt * GlobalMethods.dx * GlobalMethods.dx; }
                                    //
                                    if (guiVariables.Daily_water)
                                    {
                                        eluviated_kg *= waterfactor[row, col];
                                        ;

                                    }

                                    if (eluviated_kg > GlobalMethods.texture_kg[row, col, layer, 4]) { eluviated_kg = GlobalMethods.texture_kg[row, col, layer, 4]; }


                                    total_fine_eluviated_mass_kg += eluviated_kg;
                                    GlobalMethods.texture_kg[row, col, layer, 4] -= eluviated_kg;
                                    GlobalMethods.texture_kg[row, col, layer + 1, 4] += eluviated_kg;
                                    //improve for lowers - where does the fine clay go?
                                    //count the amount of clay and leached chem exiting catchment
                                    // SOIL possibly improve with coarse clay fraction
                                    depth += GlobalMethods.layerthickness_m[row, col, layer] / 2;
                                }
                            }
                        }
                    }
                }
                if (guiVariables.Timeseries.Total_fine_eluviated_checkbox)
                {
                    timeseries_matrix[GlobalMethods.t, timeseries_order[18]] = total_fine_eluviated_mass_kg;
                }
            }
            catch { Debug.WriteLine(" Problem occurred in translocation calculation"); }
        }

        void soil_clay_translocation_Jagercikova()
        {
            double ct_adv0, ct_adv0_all, ct_dd, ct_dd_all;
            ct_adv0_all = Convert.ToDouble(ct_v0_Jagercikova.Text);
            ct_dd_all = Convert.ToDouble(ct_dd_Jagercikova.Text);
            ct_adv0 = ct_adv0_all;
            ct_dd = ct_dd_all;

            try
            {
                // based on the advection-diffusion equation of Jagercikova et al., 2017 https://doi.org/10.1007/s11368-016-1560-9
                // We only added the advection part, because the diffusion represents bioturbation and that is already modeled elsewhere
                double local_I;


                double depth, f_clay, f_oc, d_depth, ct_advi, eluviated_kg, CEC_ct, CCEC_ct, wdclay;
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (GlobalMethods.dtm[row, col] != -9999)
                        {
                            if (NA_in_soil(row, col))
                            {
                                Debug.WriteLine("ctj1");
                            }

                            if (guiVariables.Daily_water)
                            {
                                local_I = Math.Max(Iy[row, col], 0);
                                ct_adv0 = ct_adv0_all * (1 - Math.Exp(-local_I / (2.0 * (1.0 / 3)))); // Exponential function to determine v0, based on infiltration. The function approaches a v0 of 1. I of 0.5~v0 of 0.5. the 2 indicates the range of the variogram. 
                                ct_dd = ct_dd_all - (1 - Math.Exp(-local_I / (2.0 * (1.0 / 3)))); // adjust depth decay, by subtracting 
                            }

                            depth = 0;
                            for (int layer = 0; layer < GlobalMethods.max_soil_layers; layer++) // we loop through all layers. Lowest layer has no recipient, so there we have free drainage of clay
                            {
                                if (GlobalMethods.layerthickness_m[row, col, layer] > 0)  // source layer has to exist. Adjusted for free drainage, receiving layer doesn'GlobalMethods.t have to be present
                                {

                                    if (GlobalMethods.texture_kg[row, col, layer, 3] > 0)
                                    {
                                        depth += GlobalMethods.layerthickness_m[row, col, layer] / 2;



                                        f_clay = GlobalMethods.texture_kg[row, col, layer, 3] / (GlobalMethods.texture_kg[row, col, layer, 1] + GlobalMethods.texture_kg[row, col, layer, 2] + GlobalMethods.texture_kg[row, col, layer, 3]); // fine earth fraction of clay. No fine clay
                                        f_oc = (GlobalMethods.young_SOM_kg[row, col, layer] + GlobalMethods.old_SOM_kg[row, col, layer]) / (GlobalMethods.young_SOM_kg[row, col, layer] + GlobalMethods.old_SOM_kg[row, col, layer] + GlobalMethods.texture_kg[row, col, layer, 1] + GlobalMethods.texture_kg[row, col, layer, 2] + GlobalMethods.texture_kg[row, col, layer, 3]); // fine earth fraction of clay. No fine clay
                                        f_oc /= 1.72; // calculate from SOM to SOC: https://www.researchgate.net/post/How_can_I_convert_percent_soil_organic_matter_into_soil_C

                                        if ((layer + 1) < GlobalMethods.max_soil_layers)
                                        {
                                            d_depth = (GlobalMethods.layerthickness_m[row, col, layer] + GlobalMethods.layerthickness_m[row, col, layer + 1]) / 2; // distance from mid-point to mid-point of source and sink cell
                                        }
                                        else // eluviation from lowest layer
                                        {
                                            d_depth = (GlobalMethods.layerthickness_m[row, col, layer] + GlobalMethods.layerthickness_m[row, col, layer - 1]) / 2; // use distance to higher cell as reference
                                        }

                                        // Eluviation limited by association with OM and CEC (equations from Model 2 of Brubaker et al, 1992: estimating the water-dispersible clay content of soils)
                                        // CEC estimated with PTF from Foth and Ellis 1996, as used in Finke 2012
                                        CEC_ct = (32 + 3670 * f_oc + 196 * f_clay) / 10; // cmol+/kg
                                        CCEC_ct = CEC_ct - 300 * f_oc; // carbon corrected CEC
                                        if (f_clay == 0) { f_clay = 0.000001; } // prevent dividing by 0. clay percentage of 1% always has absent dispersible clay
                                        wdclay = (0.369 * (f_clay * 100) - 8.96 * (CCEC_ct / (f_clay * 100)) + 4.48) / 100; // fraction of clay that can be dispersed
                                        if (wdclay < 0) { wdclay = 0; } // prevent negative water-dispersible clay


                                        // if (GlobalMethods.t == 3000) { Debugger.Break(); }
                                        // advection
                                        ct_advi = ct_adv0 * Math.Exp(-ct_dd * depth);
                                        eluviated_kg = ct_advi / 100 * GlobalMethods.bulkdensity[row, col, layer] * f_clay * GlobalMethods.dx * GlobalMethods.dx;
                                        // eluviated_kg = 1 / d_depth * (ct_advi * f_clay) * 1000 / GlobalMethods.bulkdensity[row, col, layer];

                                        if (eluviated_kg > (GlobalMethods.texture_kg[row, col, layer, 3] * wdclay))
                                        {
                                            eluviated_kg = GlobalMethods.texture_kg[row, col, layer, 3] * wdclay;
                                        }


                                        if (eluviated_kg > GlobalMethods.texture_kg[row, col, layer, 3]) { eluviated_kg = GlobalMethods.texture_kg[row, col, layer, 3]; } // correct for too muchy clay eluviating, not necessary anymore due to limitation water-dispersible clay




                                        total_fine_eluviated_mass_kg += eluviated_kg;
                                        GlobalMethods.texture_kg[row, col, layer, 3] -= eluviated_kg;
                                        if ((layer + 1) < GlobalMethods.max_soil_layers) { GlobalMethods.texture_kg[row, col, layer + 1, 3] += eluviated_kg; } // in case there is a lower receiving layer


                                        depth += GlobalMethods.layerthickness_m[row, col, layer] / 2;

                                    }

                                    if (NA_in_soil(row, col))
                                    {
                                        Debug.WriteLine("err_ctj2");
                                    }


                                }
                            }
                        }
                    }
                    if (guiVariables.Timeseries.Total_fine_eluviated_checkbox)
                    {
                        timeseries_matrix[GlobalMethods.t, timeseries_order[18]] = total_fine_eluviated_mass_kg;
                    }
                }
                if (NA_in_map(GlobalMethods.dtm) > 0 | NA_in_map(GlobalMethods.soildepth_m) > 0)
                {
                    Debug.WriteLine("err_ctj3");
                }

            }
            catch { Debug.WriteLine(" Problem occurred in translocation calculation"); }
        }

        void soil_silt_translocation()
        {
            //in Spitsbergen, it is mostly silt (with attendant clay) that gets translocated in the profile. Clay is not modelled in itself

            int layer;
            double eluviated_kg;
            try
            {
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        for (layer = 0; layer < GlobalMethods.max_soil_layers - 1; layer++)   // we loop through all layers except the lower one - clay translocation there has no lower recipient
                        {
                            if (GlobalMethods.layerthickness_m[row, col, layer] > 0 && GlobalMethods.layerthickness_m[row, col, layer + 1] > 0)  // both source and sink layers have to exist.
                            {
                                if (GlobalMethods.texture_kg[row, col, layer, 2] > 0)
                                {
                                    //calculate the mass of eluviation
                                    eluviated_kg = max_eluviation * (1 - Math.Exp(-Cclay * GlobalMethods.texture_kg[row, col, layer, 2])) * dt * GlobalMethods.dx * GlobalMethods.dx;
                                    GlobalMethods.texture_kg[row, col, layer, 2] -= eluviated_kg;
                                    if (GlobalMethods.texture_kg[row, col, layer, 2] < 0) { Debug.WriteLine("error: too much clay eluviation "); }
                                    GlobalMethods.texture_kg[row, col, layer + 1, 2] += eluviated_kg;
                                    //improve for lowers
                                    //count the amount of clay and leached chem exiting catchment
                                    // SOIL possibly improve with coarse clay fraction
                                }
                            }
                        }
                    }
                }
            }
            catch { Debug.WriteLine(" Problem occurred in translocation calculation"); }
        }

        void soil_decalcification()
        {
            // develop: erosion of carbonates, link to clay fraction? Or transport CO3_kg with the rest of the sediments?

            // Decalcification depends on the amount of percolation, according to Egli and Fitze (2001). The function below is a linear regression between the data in their paper. This function should work as a simple test. more complicated functions, with equilibria and secondary carbonates are possible
            try
            {
                double CO3_loss;
                // Carbonate losses [mol m-2 y-1] = 205.58 * percolation [m] - 12.392
                // Infiltration / percolation is modeled in m, so adjustments have to be made for cell size. In every step? 
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (GlobalMethods.dtm[row, col] != -9999)
                        {
                            CO3_loss = (205.58 * Iy[row, col] - 12.392) * (GlobalMethods.dx * GlobalMethods.dx) * 60.01; // Corrected the equation for cell size (GlobalMethods.dx*GlobalMethods.dx) and molar mass (60.01 g mol-1)
                            if (CO3_loss < 0) { CO3_loss = 0; }

                            int layer = 0;
                            while (CO3_loss > 0)
                            {
                                //Debug.WriteLine("dec1");
                                if (CO3_kg[row, col, layer] > 0)
                                {
                                    //Debug.WriteLine("dec1a");
                                    if (CO3_kg[row, col, layer] >= CO3_loss)
                                    {
                                        //Debug.WriteLine("dec2");
                                        CO3_kg[row, col, layer] -= CO3_loss;
                                        CO3_loss = 0;
                                    }
                                    else
                                    {
                                        //Debug.WriteLine("dec3");
                                        CO3_loss -= CO3_kg[row, col, layer];
                                        CO3_kg[row, col, layer] = 0;
                                        if (layer < (GlobalMethods.max_soil_layers - 1)) { layer++; }
                                        else { CO3_loss = 0; }// all CO3 is removed from the catchment
                                    }
                                    //Debug.WriteLine("dec4");
                                }
                                else
                                {
                                    //Debug.WriteLine("dec5");
                                    if (layer < (GlobalMethods.max_soil_layers - 1)) { layer++; }
                                    else { CO3_loss = 0; }// all CO3 is removed from the catchment
                                    ;
                                }
                            }
                            //Debug.WriteLine("layer decalc: {0}", layer);
                        }
                    }
                }
            }
            catch
            {
                MessageBox.Show("error in decalcification");
            }
        }

        void lessivage_calibration(int row, int col, int cal)
        {
            //Debug.Write(Cclay + " " + max_eluviation + " " + ct_depthdec + " ");
            //double[] lp4 = new double[,]{0.09, 0.09,0.09, 0.10, 0.16, 0.16, 0.19, 0.19, 0.19, 0.16, 0.16, 0.16, 0.16, 0.13, 0.13, 0.13, 0.13, 0.13, 0.13, 0.13};

            double[,] lp4 = new double[6, 2] { { .31, .09 }, { .45, .10 }, { .62, .16 }, { .90, .19 }, { 1.35, .16 }, { 2.0, .13 } };
            double depth, err, rmse_ct, me_ct, lp4_clay;
            int lp4_row = 0;
            int layercount = 0;
            rmse_ct = 0;
            me_ct = 0;
            depth = 0;
            Debug.WriteLine(rmse_ct + ", " + me_ct);
            for (int layer = 0; layer < GlobalMethods.max_soil_layers; layer++)
            {
                if (GlobalMethods.layerthickness_m[row, col, layer] > 0)
                {
                    depth += GlobalMethods.layerthickness_m[row, col, layer] / 2;
                    if (depth <= lp4[5, 0])
                    {
                        double totalweight = GlobalMethods.texture_kg[row, col, layer, 1] + GlobalMethods.texture_kg[row, col, layer, 2] + GlobalMethods.texture_kg[row, col, layer, 3] + GlobalMethods.texture_kg[row, col, layer, 4]; // calibrate on fine soil fraction only
                        while (depth > lp4[lp4_row, 0])
                        {
                            lp4_row++;
                        }
                        lp4_clay = lp4[lp4_row, 1];

                        err = ((GlobalMethods.texture_kg[row, col, layer, 3] + GlobalMethods.texture_kg[row, col, layer, 4]) / totalweight) - lp4_clay;
                        rmse_ct += err * err;
                        me_ct += err;
                        layercount += 1;
                    }

                    depth += GlobalMethods.layerthickness_m[row, col, layer] / 2;
                }

            }
            Debug.WriteLine(layercount);
            rmse_ct = Math.Pow(rmse_ct / layercount, .5);
            me_ct = me_ct / layercount;
            //Debug.Write(rmse_ct + " " + me_ct);
            //Debug.WriteLine("");//start on new line
            lessivage_errors[cal, 0] = Cclay;
            lessivage_errors[cal, 1] = max_eluviation;
            lessivage_errors[cal, 2] = ct_depthdec;
            lessivage_errors[cal, 3] = rmse_ct;
            lessivage_errors[cal, 4] = me_ct;
        }

        #endregion

        #region Geomorphic processes code

        int NA_in_map(double[,] map)
        {
            int NA_count = 0;

            try
            {
                for (int rowNA = 0; rowNA < GlobalMethods.nr; rowNA++)
                {
                    for (int colNA = 0; colNA < GlobalMethods.nc; colNA++)
                    {
                        if (Double.IsNaN(map[rowNA, colNA]) | Double.IsInfinity(map[rowNA, colNA]))
                        {
                            NA_count++;
                            Debug.WriteLine("NA at row {0}, col {1}", rowNA, colNA);
                        }
                    }
                }

            }
            catch
            {
                Debug.WriteLine("err_NAmap1");

            }
            return (NA_count);
        }

        bool NA_in_soil(int rowNA, int colNA)
        {
            bool boolNA = false;
            for (int layNA = 0; layNA < GlobalMethods.max_soil_layers; layNA++)
            {
                for (int texNA = 0; texNA < 5; texNA++)
                {
                    if (Double.IsNaN(GlobalMethods.texture_kg[rowNA, colNA, layNA, texNA]) | Double.IsInfinity(GlobalMethods.texture_kg[rowNA, colNA, layNA, texNA]))
                    {
                        boolNA = true;
                        Debug.WriteLine("NA in row {0}, col {1}, lay {2}, tex {3}", rowNA, colNA, layNA, texNA);
                    }

                }
                if (Double.IsNaN(GlobalMethods.young_SOM_kg[rowNA, colNA, layNA]) | Double.IsInfinity(GlobalMethods.young_SOM_kg[rowNA, colNA, layNA]))
                {
                    boolNA = true;
                    Debug.WriteLine("NA in row {0}, col {1}, lay {2}, young OM, val {3}", rowNA, colNA, layNA, GlobalMethods.young_SOM_kg[rowNA, colNA, layNA]);
                }

                if (Double.IsNaN(GlobalMethods.old_SOM_kg[rowNA, colNA, layNA]) | Double.IsInfinity(GlobalMethods.old_SOM_kg[rowNA, colNA, layNA]))
                {
                    boolNA = true;
                    Debug.WriteLine("NA in row {0}, col {1}, lay {2}, old OM, val {3}", rowNA, colNA, layNA, GlobalMethods.old_SOM_kg[rowNA, colNA, layNA]);
                }
            }

            return (boolNA);
        }

        bool NA_anywhere_in_soil()
        {

            bool boolNA = false;
            try
            {
                for (int rowNA = 0; rowNA < GlobalMethods.nr; rowNA++)
                {
                    for (int colNA = 0; colNA < GlobalMethods.nc; colNA++)
                    {
                        if (NA_in_soil(rowNA, colNA) == true) { boolNA = true; }
                    }
                }

            }
            catch
            {
                Debug.WriteLine("err_NAmap1");

            }
            return (boolNA);
        }

        int NA_in_location(double[,] map, int rowNA, int colNA)
        {
            int NA_count = 0;
            try
            {
                if (Double.IsNaN(map[rowNA, colNA]) | Double.IsInfinity(map[rowNA, colNA]))
                {
                    NA_count++;
                    Debug.WriteLine("NA at row {0}, col {1}", rowNA, colNA);
                }
            }
            catch
            {
                Debug.WriteLine("err_nal1");

            }
            return (NA_count);
        }


        void calculate_water_ero_sed_daily()
        {
            if (NA_in_map(GlobalMethods.dtm) > 0 | NA_in_map(GlobalMethods.soildepth_m) > 0)
            {
                Debug.WriteLine("we1");
            }

            double mass_before = total_catchment_mass(), mass_after, mass_export = 0;
            //Debug.WriteLine("WE1");
            int size, dir;
            double water_out, flow_between_cells_m3_per_m, total_sediment_in_transport_kg, rock_fraction, mass_to_be_eroded, selectivity_fraction, total_ero = 0, total_dep = 0, potential_transported_amount_kg, vegetation_cover_fraction;
            double[] total_mass_eroded = new double[7] { 0, 0, 0, 0, 0, 0, 0 };
            double[] total_mass_deposited = new double[7] { 0, 0, 0, 0, 0, 0, 0 };


            double[,,] local_mass_balance = new double[GlobalMethods.nr, GlobalMethods.nc, 5]; // 0: mass in, 1: mass out, 2: mass ero, 3: mass sed, 4: Q_out

            double[,] mass_difference_input_output = new double[GlobalMethods.nr, GlobalMethods.nc];

            // 1: set all water and sediment flow to 0
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    {
                        if (only_waterflow_checkbox.Checked == false)
                        {

                            for (size = 0; size < GlobalMethods.n_texture_classes; size++)
                            {
                                GlobalMethods.sediment_in_transport_kg[row, col, size] = 0;
                            }
                            GlobalMethods.old_SOM_in_transport_kg[row, col] = 0;
                            GlobalMethods.young_SOM_in_transport_kg[row, col] = 0;
                            GlobalMethods.dz_ero_m[row, col] = 0;
                            GlobalMethods.dz_sed_m[row, col] = 0;
                            GlobalMethods.lake_sed_m[row, col] = 0;

                            for (int iter = 0; iter < 5; iter++)
                            {
                                local_mass_balance[row, col, iter] = 0;
                            }
                        }
                    }
                }  // end for col
            }  //end for row
               //Debug.WriteLine("WE2");

            // 2: Iterate through rows and columns
            int runner = 0;
            for (runner = GlobalMethods.number_of_data_cells - 1; runner >= 0; runner--)
            {
                int row, col;
                if (GlobalMethods.index[runner] != -9999)
                {
                    //if (row == 40 & col == 31) { displaysoil(40, 31); Debugger.Break(); }

                    row = GlobalMethods.row_index[runner]; col = GlobalMethods.col_index[runner];

                    // local mass balance
                    local_mass_balance[row, col, 4] = guiVariables.OFy_m[row, col, 0];

                    //Debug.WriteLine("WE3");

                    // 3: Determine fraction of water flowing to lower neighbour
                    water_out = 0;
                    for (int dir2 = 1; dir2 < 9; dir2++)
                    {
                        water_out += guiVariables.OFy_m[row, col, dir2];
                        //if (col == 12 && dir2 == 1) { Debugger.Break(); }
                        // if (GlobalMethods.t == 1) { Debugger.Break(); }
                    }
                    double fracsum = 0;


                    if (water_out > 0)
                    {
                        // Debug.WriteLine("Overland flow in col {0}", col);
                        dir = 0;
                        for (i = (-1); i <= 1; i++)
                        {
                            for (j = (-1); j <= 1; j++)
                            {
                                if (!((i == 0) && (j == 0))) { dir++; }

                                dh = 0; fraction = 0; transport_capacity_kg = 0;


                                if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                                {
                                    if (GlobalMethods.dtm[row + i, col + j] != -9999)
                                    {
                                        // Debug.WriteLine("row = {0}, col = {1}, dir = {2}, i = {3}, j = {4}",row,col,dir,i,j);
                                        dh = GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + i, col + j];
                                        //Debug.WriteLine("Source: {0}, sink: {1}", GlobalMethods.dtm[row, col], GlobalMethods.dtm[row + i, col + j]);
                                        GlobalMethods.d_x = GlobalMethods.dx;
                                        //if (col + j == 0) { Debugger.Break(); }

                                        if (dh > 0)
                                        {
                                            if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                                            dh /= GlobalMethods.d_x;

                                            fraction = guiVariables.OFy_m[row, col, dir] / water_out;
                                            fracsum += fraction;
                                            flow_between_cells_m3_per_m = guiVariables.OFy_m[row, col, dir] * GlobalMethods.dx * GlobalMethods.dx / GlobalMethods.dx; // 

                                            //Debug.WriteLine("WE4");

                                            // 4: Determine transport capacity
                                            transport_capacity_kg = advection_erodibility * (GlobalMethods.bulkdensity[row, col, 0] * GlobalMethods.dx * GlobalMethods.dx) * (Math.Pow(flow_between_cells_m3_per_m, m) * Math.Pow(dh, n));
                                            if (transport_capacity_kg < 0) { transport_capacity_kg = 0; Debug.WriteLine(" Warning: negative transport capacity at" + row + " " + col); } // warning. This should never happen

                                            total_sediment_in_transport_kg = 0;
                                            for (size = 0; size < GlobalMethods.n_texture_classes; size++)
                                            {
                                                total_sediment_in_transport_kg += fraction * GlobalMethods.sediment_in_transport_kg[row, col, size]; //all in kg.
                                            } // organic matter doesn'GlobalMethods.t count for total matter in transport


                                            // 5: Determine transport / erosion / sedimentation

                                            //Debug.WriteLine("WE5a");

                                            // 5a: Transport
                                            if (transport_capacity_kg == total_sediment_in_transport_kg)
                                            {
                                                // everything gets transported
                                                for (size = 0; size < GlobalMethods.n_texture_classes; size++)
                                                {
                                                    GlobalMethods.sediment_in_transport_kg[row + i, col + j, size] += fraction * GlobalMethods.sediment_in_transport_kg[row, col, size];  //all in kg 

                                                    // local mass balance
                                                    local_mass_balance[row + i, col + j, 0] += fraction * GlobalMethods.sediment_in_transport_kg[row, col, size]; // incoming mass
                                                    local_mass_balance[row, col, 1] += fraction * GlobalMethods.sediment_in_transport_kg[row, col, size]; // outgoing mass
                                                }
                                                GlobalMethods.old_SOM_in_transport_kg[row + i, col + j] += fraction * GlobalMethods.old_SOM_in_transport_kg[row, col];  //all in kg
                                                GlobalMethods.young_SOM_in_transport_kg[row + i, col + j] += fraction * GlobalMethods.young_SOM_in_transport_kg[row, col];  //all in kg
                                            }

                                            //Debug.WriteLine("WE5b");

                                            // 5b: Erosion
                                            if (transport_capacity_kg > total_sediment_in_transport_kg)
                                            {
                                                // everything we want to transport gets transported, plus a little bit extra
                                                for (size = 0; size < GlobalMethods.n_texture_classes; size++)
                                                {
                                                    GlobalMethods.sediment_in_transport_kg[row + i, col + j, size] += fraction * GlobalMethods.sediment_in_transport_kg[row, col, size];  //all in kg 

                                                    // local mass balance
                                                    local_mass_balance[row + i, col + j, 0] += fraction * GlobalMethods.sediment_in_transport_kg[row, col, size]; // incoming mass
                                                    local_mass_balance[row, col, 1] += fraction * GlobalMethods.sediment_in_transport_kg[row, col, size]; // outgoing mass
                                                }
                                                GlobalMethods.old_SOM_in_transport_kg[row + i, col + j] += fraction * GlobalMethods.old_SOM_in_transport_kg[row, col];  //all in kg
                                                GlobalMethods.young_SOM_in_transport_kg[row + i, col + j] += fraction * GlobalMethods.young_SOM_in_transport_kg[row, col];  //all in kg

                                                // to calculate the (possibly) extra amount we are going to transport, we first evaluate whether we exceed the erosion threshold   
                                                if ((transport_capacity_kg - total_sediment_in_transport_kg) > erosion_threshold_kg)
                                                {

                                                    //first, calculate how much we are going to erode. Not as much as we want to if the soil is protected by rocks or plants
                                                    rock_fraction = GlobalMethods.texture_kg[row, col, 0, 0] / (GlobalMethods.texture_kg[row, col, 0, 0] + GlobalMethods.texture_kg[row, col, 0, 1] + GlobalMethods.texture_kg[row, col, 0, 2] + GlobalMethods.texture_kg[row, col, 0, 3] + GlobalMethods.texture_kg[row, col, 0, 4]);
                                                    if (aridity_vegetation[row, col] >= 1) { vegetation_cover_fraction = 1; }
                                                    else { vegetation_cover_fraction = aridity_vegetation[row, col]; }// 
                                                    mass_to_be_eroded = (transport_capacity_kg - total_sediment_in_transport_kg) * Math.Exp(-rock_protection_constant * rock_fraction) * Math.Exp(-bio_protection_constant * vegetation_cover_fraction);
                                                    // DEVELOP MM changing vegetation cover in a year

                                                    // second, calculate how the mass to be eroded is taken from the different size fractions: selectivity

                                                    // if total transport capacity is small, only the finer fractions will be eroded (selectivity with diameter to power 0.5). 
                                                    //For larger transport capacities, selectivity decreases (diameter to power 0 = equal between fractions)

                                                    double constant_b1 = 0.5 * Math.Exp(constant_selective_transcap * transport_capacity_kg);
                                                    double sum_diameter_power = 0;
                                                    for (size = 0; size < 5; size++)
                                                    {
                                                        sum_diameter_power += 1 / Math.Pow(upper_particle_size[size], constant_b1);
                                                    }
                                                    double clayeroded_0_kg = 0, claypresent_0_kg = 0, clayeroded_1_kg = 0, claypresent_1_kg = 0;


                                                    for (size = 0; size < 5; size++)
                                                    {

                                                        selectivity_fraction = (1 / Math.Pow(upper_particle_size[size], constant_b1)) / sum_diameter_power;    // unit [-]
                                                        if (GlobalMethods.texture_kg[row, col, 0, size] >= selectivity_fraction * mass_to_be_eroded)
                                                        {    // typical situation
                                                            if (size > 2)
                                                            {
                                                                clayeroded_0_kg += selectivity_fraction * mass_to_be_eroded;
                                                                claypresent_0_kg += GlobalMethods.texture_kg[row, col, 0, size];

                                                            }
                                                            mass_difference_input_output[row, col] += selectivity_fraction * mass_to_be_eroded;
                                                            total_mass_eroded[size] += selectivity_fraction * mass_to_be_eroded;
                                                            GlobalMethods.texture_kg[row, col, 0, size] -= selectivity_fraction * mass_to_be_eroded;   // unit [kg]
                                                            GlobalMethods.sediment_in_transport_kg[row + i, col + j, size] += selectivity_fraction * mass_to_be_eroded;  // unit [kg]
                                                            total_ero += selectivity_fraction * mass_to_be_eroded;


                                                            // local mass balance
                                                            local_mass_balance[row + i, col + j, 0] += selectivity_fraction * mass_to_be_eroded; // incoming mass next cell
                                                            local_mass_balance[row, col, 1] += selectivity_fraction * mass_to_be_eroded; // outgoing mass current cell
                                                            local_mass_balance[row, col, 2] += selectivity_fraction * mass_to_be_eroded; // eroded mass current cell


                                                        }
                                                        else
                                                        {
                                                            // exceptional. If we want to erode more than present in the layer, we will take it from one layer down.
                                                            //this is to avoid exceptionally thin rocky layers blocking all erosion
                                                            mass_difference_input_output[row, col] += GlobalMethods.texture_kg[row, col, 0, size];
                                                            total_mass_eroded[size] += GlobalMethods.texture_kg[row, col, 0, size];
                                                            double left = (selectivity_fraction * mass_to_be_eroded) - GlobalMethods.texture_kg[row, col, 0, size]; // unit [kg]

                                                            GlobalMethods.sediment_in_transport_kg[row + i, col + j, size] += GlobalMethods.texture_kg[row, col, 0, size];
                                                            total_ero += GlobalMethods.texture_kg[row, col, 0, size];

                                                            if (size > 2)
                                                            {
                                                                clayeroded_0_kg += GlobalMethods.texture_kg[row, col, 0, size];
                                                                claypresent_0_kg += 0;
                                                            }

                                                            // local mass balance
                                                            local_mass_balance[row + i, col + j, 0] += GlobalMethods.texture_kg[row, col, 0, size]; // incoming mass next cell
                                                            local_mass_balance[row, col, 1] += GlobalMethods.texture_kg[row, col, 0, size]; // outgoing mass current cell
                                                            local_mass_balance[row, col, 2] += GlobalMethods.texture_kg[row, col, 0, size]; // eroded mass current cell

                                                            GlobalMethods.texture_kg[row, col, 0, size] = 0;

                                                            if (GlobalMethods.texture_kg[row, col, 1, size] >= left)
                                                            {   // typical
                                                                mass_difference_input_output[row, col] += left;
                                                                total_mass_eroded[size] += left;
                                                                if (size > 2)
                                                                {
                                                                    clayeroded_1_kg += left;
                                                                    claypresent_1_kg += GlobalMethods.texture_kg[row, col, 1, size] - left;
                                                                }

                                                                GlobalMethods.texture_kg[row, col, 1, size] -= left;  // unit [kg]
                                                                GlobalMethods.sediment_in_transport_kg[row + i, col + j, size] += left;  // unit [kg]
                                                                total_ero += left;


                                                                // local mass balance
                                                                local_mass_balance[row + i, col + j, 0] += left; // incoming mass next cell
                                                                local_mass_balance[row, col, 1] += left; // outgoing mass current cell
                                                                local_mass_balance[row, col, 2] += left; // eroded mass current cell


                                                            }
                                                            else
                                                            {

                                                                total_mass_eroded[size] += GlobalMethods.texture_kg[row, col, 1, size];
                                                                mass_difference_input_output[row, col] += GlobalMethods.texture_kg[row, col, 1, size];

                                                                if (size > 2)
                                                                {
                                                                    clayeroded_1_kg += GlobalMethods.texture_kg[row, col, 1, size];
                                                                    claypresent_1_kg += 0; //MM Leidt dit niet tot delen door nul?
                                                                }

                                                                GlobalMethods.sediment_in_transport_kg[row + i, col + j, size] += GlobalMethods.texture_kg[row, col, 1, size];// unit [kg] // MM adjusted in water erosion
                                                                total_ero += GlobalMethods.texture_kg[row, col, 1, size];

                                                                // local mass balance
                                                                local_mass_balance[row + i, col + j, 0] += GlobalMethods.texture_kg[row, col, 1, size]; // incoming mass next cell
                                                                local_mass_balance[row, col, 1] += GlobalMethods.texture_kg[row, col, 1, size]; // outgoing mass current cell
                                                                local_mass_balance[row, col, 2] += GlobalMethods.texture_kg[row, col, 1, size]; // eroded mass current cell


                                                                GlobalMethods.texture_kg[row, col, 1, size] = 0;

                                                            } // end else 

                                                        }
                                                    } // end size

                                                    //organic matter. // 3c2. Organic matter maximaal erosie = die fractie van OM in de geerodeerde lagen die overeenkomt met de fractie geerodeerde fine earth in de geerodeerde lagen. is eroded as a fraction of total OM. That fraction equals the fraction of clay eroded from the layer
                                                    // do not forget (after all erosion and depositionthe assumption underlying this is done) to recalculate the thickness of layers and adapt the GlobalMethods.dtm to that. clay and humus are bound in aggregates
                                                    //this does not cover: LMW SOM, peat or large woody debris
                                                    double clayerodedfraction_0;
                                                    if (clayeroded_0_kg > 0) { clayerodedfraction_0 = clayeroded_0_kg / (clayeroded_0_kg + claypresent_0_kg); }
                                                    else { clayerodedfraction_0 = 0; }


                                                    double clayerodedfraction_1;
                                                    if (clayeroded_1_kg > 0) { clayerodedfraction_1 = clayeroded_1_kg / (clayeroded_1_kg + claypresent_1_kg); }
                                                    else { clayerodedfraction_1 = 0; }
                                                    GlobalMethods.old_SOM_in_transport_kg[row + i, col + j] += GlobalMethods.old_SOM_kg[row, col, 0] * clayerodedfraction_0 + GlobalMethods.old_SOM_kg[row, col, 1] * clayerodedfraction_1;
                                                    GlobalMethods.young_SOM_in_transport_kg[row + i, col + j] += GlobalMethods.young_SOM_kg[row, col, 0] * clayerodedfraction_0 + GlobalMethods.young_SOM_kg[row, col, 1] * clayerodedfraction_1;
                                                    total_mass_eroded[5] += GlobalMethods.old_SOM_kg[row, col, 0] * clayerodedfraction_0 + GlobalMethods.old_SOM_kg[row, col, 1] * clayerodedfraction_1;
                                                    total_mass_eroded[6] += GlobalMethods.young_SOM_kg[row, col, 0] * clayerodedfraction_0 + GlobalMethods.young_SOM_kg[row, col, 1] * clayerodedfraction_1;
                                                    GlobalMethods.old_SOM_kg[row, col, 0] *= (1 - clayerodedfraction_0);
                                                    GlobalMethods.young_SOM_kg[row, col, 0] *= (1 - clayerodedfraction_0);
                                                    GlobalMethods.old_SOM_kg[row, col, 1] *= (1 - clayerodedfraction_1);
                                                    GlobalMethods.young_SOM_kg[row, col, 1] *= (1 - clayerodedfraction_1);

                                                    if (double.IsNaN(GlobalMethods.young_SOM_kg[row, col, 0]))
                                                    {
                                                        Debug.WriteLine("err_we2");
                                                    }

                                                    //Local Mass balance klopt niet meer met OM erbij



                                                }
                                                else
                                                {
                                                    //do nothing. We wanted to erode, but not enough so to actually exceed the threshold and actually do that
                                                }
                                            }
                                            //Debug.WriteLine("WE5c");

                                            // 5c: Sedimentation
                                            if (transport_capacity_kg < total_sediment_in_transport_kg)
                                            {

                                                //first, calculate how much we are going to keep in transport. This is the way that selectivity works now. 
                                                double sum_diameter_power = 0, clay_deposited = 0, clay_transported = 0;
                                                for (size = 0; size < 5; size++)
                                                {
                                                    sum_diameter_power += 1 / Math.Pow(upper_particle_size[size], 0.5);
                                                }

                                                for (size = 0; size < 5; size++)
                                                {
                                                    selectivity_fraction = (1 / Math.Pow(upper_particle_size[size], 0.5)) / sum_diameter_power;    // unit [-]
                                                    potential_transported_amount_kg = selectivity_fraction * transport_capacity_kg;                      // unit [kg]

                                                    if (potential_transported_amount_kg < (fraction * GlobalMethods.sediment_in_transport_kg[row, col, size]))
                                                    {
                                                        total_mass_deposited[size] += fraction * GlobalMethods.sediment_in_transport_kg[row, col, size] - potential_transported_amount_kg;
                                                        total_mass_eroded[size] -= fraction * GlobalMethods.sediment_in_transport_kg[row, col, size] - potential_transported_amount_kg;
                                                        mass_difference_input_output[row, col] -= fraction * GlobalMethods.sediment_in_transport_kg[row, col, size] - potential_transported_amount_kg;

                                                        // local mass balance
                                                        local_mass_balance[row, col, 3] += fraction * GlobalMethods.sediment_in_transport_kg[row, col, size] - potential_transported_amount_kg; // deposition mass current cell
                                                        local_mass_balance[row + i, col + j, 0] += potential_transported_amount_kg; // incoming mass next cell
                                                        local_mass_balance[row, col, 1] += potential_transported_amount_kg; // outgoing mass current cell


                                                        GlobalMethods.texture_kg[row, col, 0, size] += fraction * GlobalMethods.sediment_in_transport_kg[row, col, size] - potential_transported_amount_kg;        // unit [kg]
                                                        total_dep += fraction * GlobalMethods.sediment_in_transport_kg[row, col, size] - potential_transported_amount_kg;
                                                        GlobalMethods.sediment_in_transport_kg[row + i, col + j, size] += potential_transported_amount_kg;                                    // unit [kg]  


                                                        if (size > 2)
                                                        {
                                                            clay_deposited += fraction * GlobalMethods.sediment_in_transport_kg[row, col, size] - potential_transported_amount_kg;
                                                            clay_transported += potential_transported_amount_kg;
                                                        }

                                                        //if(size==3 || size==4) { organic_selectivity_fraction += selectivity_fraction; }
                                                    }
                                                    else
                                                    {
                                                        //do nothing. We keep the sediment in transport, and do not deposit anything
                                                        GlobalMethods.sediment_in_transport_kg[row + i, col + j, size] += fraction * GlobalMethods.sediment_in_transport_kg[row, col, size];


                                                        // local mass balance
                                                        local_mass_balance[row + i, col + j, 0] += fraction * GlobalMethods.sediment_in_transport_kg[row, col, size]; // incoming mass next cell
                                                        local_mass_balance[row, col, 1] += fraction * GlobalMethods.sediment_in_transport_kg[row, col, size]; // outgoing mass current cell
                                                    }
                                                }
                                                // now organic matter. 
                                                double claydepfraction = 0;
                                                if (clay_deposited > 0) { claydepfraction = clay_deposited / (clay_deposited + clay_transported); }
                                                total_mass_deposited[5] += fraction * GlobalMethods.young_SOM_in_transport_kg[row, col] * claydepfraction;
                                                total_mass_deposited[6] += fraction * GlobalMethods.old_SOM_in_transport_kg[row, col] * claydepfraction;
                                                GlobalMethods.young_SOM_kg[row, col, 0] += fraction * GlobalMethods.young_SOM_in_transport_kg[row, col] * claydepfraction;
                                                GlobalMethods.old_SOM_kg[row, col, 0] += fraction * GlobalMethods.old_SOM_in_transport_kg[row, col] * claydepfraction;

                                                GlobalMethods.young_SOM_in_transport_kg[row + i, col + j] += fraction * GlobalMethods.young_SOM_in_transport_kg[row, col] * (1 - claydepfraction);
                                                GlobalMethods.old_SOM_in_transport_kg[row + i, col + j] += fraction * GlobalMethods.old_SOM_in_transport_kg[row, col] * (1 - claydepfraction);
                                                //MM aangepast, om ook OM naar de volgende cel te transporteren 


                                                //organic_selectivity_fraction /= 2;  // unit [-]
                                                // the above procedure may cause instability problems by depositing lots of material
                                            } // end sedimentation
                                        } // end dh > 0
                                    } // end GlobalMethods.dtm!=-9999
                                }
                            } // end j
                        } // end i

                        if (fracsum < 0.9999 & !search_nodataneighbour(row, col))
                        {
                            Debug.WriteLine("fracsum = " + fracsum);
                            for (int otp = 0; otp < 10; otp++)
                            {
                                Debug.WriteLine("dir {0}, {1}", otp, guiVariables.OFy_m[row, col, otp]);
                            }
                            //minimaps(row, col);
                            Debug.WriteLine("err_we3");

                        }
                        // GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + GlobalMethods.t + "_out_dtm_test.asc", GlobalMethods.dtm);

                    } // outflow water
                    else
                    { // no outflow of water. 
                      // -> Cell located in depression. Deposit all sediments
                      // -> Cell at border of landscape, outflow of all sediments
                      // double mass_temp = total_catchment_mass();
                        bool bool_outflow = false;
                        for (i = (-1); i <= 1; i++)
                        {
                            for (j = (-1); j <= 1; j++)
                            {
                                if (i != 0 & j != 0)
                                {
                                    if ((row + i) >= GlobalMethods.nr | (row + i) <= 0 | (col + j) >= GlobalMethods.nc | (col + j) <= 0) // Does the cell fall outside of the area?
                                    {
                                        bool_outflow = true;
                                    }
                                    else if (GlobalMethods.dtm[row + i, col + j] == -9999)
                                    {
                                        bool_outflow = true;
                                    }
                                }
                            }
                        }
                        if (bool_outflow) // if there is outflow, export of sediments. 
                        {
                            for (size = 0; size < GlobalMethods.n_texture_classes; size++)
                            {
                                mass_export += GlobalMethods.sediment_in_transport_kg[row, col, size];
                            }
                            mass_export += GlobalMethods.old_SOM_in_transport_kg[row, col];  //all in kg
                            mass_export += GlobalMethods.young_SOM_in_transport_kg[row, col];  //all in kg
                        }
                        else // if there is no outflow, deposition of sediments in cell. No delta formation (yet) develop
                        {
                            for (size = 0; size < GlobalMethods.n_texture_classes; size++)
                            {

                                total_mass_deposited[size] += GlobalMethods.sediment_in_transport_kg[row, col, size];
                                total_mass_eroded[size] -= GlobalMethods.sediment_in_transport_kg[row, col, size];
                                mass_difference_input_output[row, col] -= GlobalMethods.sediment_in_transport_kg[row, col, size];

                                GlobalMethods.texture_kg[row, col, 0, size] += GlobalMethods.sediment_in_transport_kg[row, col, size];  //all in kg 
                                total_dep += GlobalMethods.sediment_in_transport_kg[row, col, size];

                                // local mass balance
                                local_mass_balance[row, col, 3] += GlobalMethods.sediment_in_transport_kg[row, col, size];
                            }
                            GlobalMethods.old_SOM_kg[row, col, 0] += GlobalMethods.old_SOM_in_transport_kg[row, col];  //all in kg
                            GlobalMethods.young_SOM_kg[row, col, 0] += GlobalMethods.young_SOM_in_transport_kg[row, col];  //all in kg
                        }


                        //displaysoil(row, col); Debugger.Break();


                        //if (row == 40 & col == 31) { displaysoil(40, 31); Debugger.Break(); }
                        // double mass_temp_diff = mass_temp - total_catchment_mass();
                    }
                    //Debug.WriteLine("WE6");

                    double mass_bal = (local_mass_balance[row, col, 0] + local_mass_balance[row, col, 2]) - (local_mass_balance[row, col, 1] + local_mass_balance[row, col, 3]); // (transport_in + depo) - (transport_out + ero)
                                                                                                                                                                                 // Debug.WriteLine("GlobalMethods.t: {0}, row: {1}, col: {2}, mass balance {3}",GlobalMethods.t, row, col,mass_bal);
                                                                                                                                                                                 // Debug.WriteLine("GlobalMethods.t: {0}, row: {1}, col: {2}, t_in: {3}, t_out: {4}, ero: {5}, depo: {6}, Q: {7}, pedonmass: {8}", GlobalMethods.t, row, col, local_mass_balance[row, col, 0], local_mass_balance[row, col, 1], local_mass_balance[row, col, 2], local_mass_balance[row, col, 3], local_mass_balance[row, col, 4],total_soil_mass(row, col));
                }
            }
            mass_after = total_catchment_mass();
            if (mass_before - (mass_after + mass_export) > 0.001)
            {
                Debug.WriteLine("err_we4");
            }
            // if (Math.Round(mass_before,6) != Math.Round(mass_after + mass_export,6)) { Debugger.Break(); }


            // 6: Calculate elevation change
            // all cells have now been considered in order of (original) altitude. We must still recalculate their thicknesses and recalculate altitude. While doing that, we should count how much erosion and deposition there has been. 
            double old_total_elevation = total_catchment_elevation();
            volume_eroded = 0; sediment_exported = 0; volume_deposited = 0;
            total_average_altitude = 0; total_altitude = 0;
            total_rain = 0; total_evap = 0; total_infil = 0; total_outflow = 0;
            wet_cells = 0; eroded_cells = 0; deposited_cells = 0;
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    if (GlobalMethods.dtm[row, col] != -9999)
                    {
                        if (only_waterflow_checkbox.Checked == false)
                        {
                            double old_thickness = GlobalMethods.soildepth_m[row, col];
                            update_all_soil_thicknesses(row, col);
                            double new_thickness = total_soil_thickness(row, col);

                            GlobalMethods.dtm[row, col] += new_thickness - old_thickness;
                            GlobalMethods.soildepth_m[row, col] = new_thickness;

                            GlobalMethods.dtmchange[row, col] += new_thickness - old_thickness;
                            GlobalMethods.sum_water_erosion[row, col] += new_thickness - old_thickness;

                            if ((new_thickness - old_thickness) < 0)
                            { GlobalMethods.dz_ero_m[row, col] += new_thickness - old_thickness; }
                            else { GlobalMethods.dz_sed_m[row, col] += new_thickness - old_thickness; }

                            if (-GlobalMethods.dz_ero_m[row, col] > guiVariables.Timeseries.Timeseries_erosion_threshold) { eroded_cells++; }
                            if (GlobalMethods.dz_sed_m[row, col] + GlobalMethods.lake_sed_m[row, col] > guiVariables.Timeseries.Timeseries_deposition_threshold) { deposited_cells++; }
                        }

                        // 7: Update timeseries
                        if (check_space_rain.Checked == true) { total_rain += GlobalMethods.rain[row, col]; }
                        total_rain += rain_value_m;
                        if (check_space_evap.Checked == true) { total_evap += GlobalMethods.evapotranspiration[row, col]; }
                        total_evap += evap_value_m;
                        if (check_space_infil.Checked == true) { total_infil += GlobalMethods.infil[row, col]; }
                        total_infil += infil_value_m;
                        if (GlobalMethods.waterflow_m3[row, col] * GlobalMethods.dx * GlobalMethods.dx > guiVariables.Timeseries.Timeseries_waterflow_threshold) { wet_cells++; }
                    } // end for nodata
                }   // end for col
            } // end for row

            //double new_total_elevation = total_catchment_elevation();
            //if(Math.Abs(new_total_elevation-old_total_elevation)>0.001)
            //{
            //    int t_err = GlobalMethods.t;
            //    Debugger.Break();
            //}


            // GlobalMethods.out_double(GlobalMethods.Workdir + "\\" + run_number + "_" + GlobalMethods.t + "_mass_difference.asc", mass_difference_input_output);
            total_rain *= GlobalMethods.dx * GlobalMethods.dx;   // m3
            total_evap *= GlobalMethods.dx * GlobalMethods.dx;   // m3
            total_infil *= GlobalMethods.dx * GlobalMethods.dx;  // m3
            total_outflow = total_rain - total_evap - total_infil;
            //Debug.WriteLine("\n--erosion and deposition overview--");
            //Debug.WriteLine("GlobalMethods.rain " + total_rain + " evap " + total_evap + " total_infil " + total_infil);
            if (only_waterflow_checkbox.Checked == false)
            {
                double total_kg_eroded = 0, total_kg_deposited = 0;
                for (size = 0; size < 7; size++)
                {
                    total_kg_eroded += total_mass_eroded[size];
                    total_kg_deposited += total_mass_deposited[size];
                }
            }
            this.InfoStatusPanel.Text = "calc movement has been finished";
            this.out_sed_statuspanel.Text = string.Format("sed_exp {0:F0} * 1000 m3", total_sed_export * GlobalMethods.dx * GlobalMethods.dx / 1000);

            //double mass_after = total_catchment_mass();
            //if(mass_before != mass_after)
            //{
            //    Debug.WriteLine("");
            //    Debug.Write(mass_before - mass_after);
            //    Debug.Write(",");
            //    Debug.Write(total_ero - total_dep);
            //    Debug.Write(",");
            //    Debug.Write(total_ero);
            //    Debug.Write(",");
            //    Debug.Write(total_mass_eroded.Sum());
            //    Debug.Write(",");

            //    // Debug.Write(total_sediment_in_transport_kg());
            //    // Debugger.Break(); 
            //}
            if (NA_in_map(GlobalMethods.dtm) > 0 | NA_in_map(GlobalMethods.soildepth_m) > 0)
            {
                Debug.WriteLine("err_we5");
            }

        }


        void calculate_water_ero_sed()    //where the water starts flowing, eroding and transporting
        {
            this.InfoStatusPanel.Text = "water erosion calculation";
            dhmax_errors = 0;
            //set all start q values effective precipitation at time GlobalMethods.t
            nb_ok = 0;  // nb_ok is 1 als er uberhaupt buren zijn, dus 0 als er alleen maar NODATA is
            nb_check = 0; all_grids = 0;
            dz_bal = 0; sediment_exported = 0; erocnt = 0; sedcnt = 0;
            sedbal = 0; erobal = 0; maximum_allowed_deposition = -9999.0; dh_tol = 0.00025;
            sedbal2 = 0; erobal2 = 0;
            tel1 = 0; tel2 = 0; tel3 = 0; tel4 = 0;
            depressions_filled = 0; depressions_delta = 0; depressions_alone = 0; sediment_delta_m = 0; sediment_filled_m = 0; depressionvolume_filled_m = 0; crashed = false;

            double powered_slope_sum, flow_between_cells_m3_per_m, total_sediment_in_transport_kg, organic_in_transport, mass_to_be_eroded, rock_fraction, bio_fraction, selectivity_fraction, potential_transported_amount_kg, organic_selectivity_fraction;
            int size;
            double[] total_mass_eroded, total_mass_deposited_kg;
            total_mass_eroded = new double[7] { 0, 0, 0, 0, 0, 0, 0 };
            total_mass_deposited_kg = new double[7] { 0, 0, 0, 0, 0, 0, 0 };

            for (alpha = 1; alpha <= maxdepressionnumber; alpha++)  // zeroing all waterflow at outlets of depressions
            {
                depressionconsidered[alpha] = 0;
                for (int outletcounter = 0; outletcounter < 5; outletcounter++)
                {
                    if (GlobalMethods.drainingoutlet_col[alpha, outletcounter] != -1)
                    {
                        GlobalMethods.waterflow_m3[GlobalMethods.drainingoutlet_col[alpha, outletcounter], GlobalMethods.drainingoutlet_col[alpha, outletcounter]] = 0;
                    }
                }
            }
            if (NA_anywhere_in_soil() == true) { Debug.WriteLine("NA found before row col loop in water erosed"); }
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {

                    if (GlobalMethods.dtm[row, col] != -9999)
                    {
                        // First, we apply rainwater to our landscape (in a two step approach - first normal cells and lake outlets)
                        if (GlobalMethods.depression[row, col] == 0 ||
                            (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 0] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 0] == col) ||
                            (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 1] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 1] == col) ||
                            (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 2] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 2] == col) ||
                            (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 3] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 3] == col) ||
                            (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 4] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 4] == col))
                        {
                            if (check_space_evap.Checked == true) { evap_value_m = GlobalMethods.evapotranspiration[row, col]; }
                            if (check_space_rain.Checked == true) { rain_value_m = GlobalMethods.rain[row, col]; }
                            if (check_space_infil.Checked == true) { infil_value_m = GlobalMethods.infil[row, col]; }
                            //ArT // development required to account for f(GlobalMethods.t) situations
                            GlobalMethods.waterflow_m3[row, col] += (rain_value_m - infil_value_m - evap_value_m) * GlobalMethods.dx * GlobalMethods.dx;
                            if (GlobalMethods.waterflow_m3[row, col] < 0) { GlobalMethods.waterflow_m3[row, col] = 0; }
                            if (GlobalMethods.waterflow_m3[row, col] < -0.001) { Debug.WriteLine(" Negative waterflow at " + row + " " + col + ": " + GlobalMethods.waterflow_m3[row, col] + ". GlobalMethods.rain " + rain_value_m + " GlobalMethods.infil " + infil_value_m + " evap " + evap_value_m + " use " + GlobalMethods.landuse[row, col]); }
                        }
                        else  // and then, second step, for other lake cells
                        { // for other lakecells, we send the rainwater directly (equally distributed) to that lake's outlet(s) (infiltration is not zero in the lake at the moment)
                            //Debug.WriteLine(" B at " + row + " col " + col + " alt " + GlobalMethods.dtm[row, col] + " dep " + depression[row, col]);
                            int outletcounter = 0; ;
                            while (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], outletcounter] != -1)
                            {
                                outletcounter++;
                                if (outletcounter == 5) { break; }
                            }
                            for (i = 0; i < outletcounter; i++)
                            {

                                if (check_space_evap.Checked == true) { evap_value_m = GlobalMethods.evapotranspiration[row, col]; }
                                if (check_space_rain.Checked == true) { rain_value_m = GlobalMethods.rain[row, col]; }
                                if (check_space_infil.Checked == true) { infil_value_m = GlobalMethods.infil[row, col]; }
                                //ArT // development required to account for f(GlobalMethods.t) situations
                                //ArT remember to check for negative lake outflow once it happens
                                GlobalMethods.waterflow_m3[GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], i], GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], i]] += GlobalMethods.dx * GlobalMethods.dx * (rain_value_m - infil_value_m - evap_value_m) / outletcounter;
                            }
                        }
                        if (only_waterflow_checkbox.Checked == false)
                        {
                            for (size = 0; size < GlobalMethods.n_texture_classes; size++)
                            {
                                GlobalMethods.sediment_in_transport_kg[row, col, size] = 0;
                            }
                            GlobalMethods.old_SOM_in_transport_kg[row, col] = 0;
                            GlobalMethods.young_SOM_in_transport_kg[row, col] = 0;
                            GlobalMethods.dz_ero_m[row, col] = 0;
                            GlobalMethods.dz_sed_m[row, col] = 0;
                            GlobalMethods.lake_sed_m[row, col] = 0;

                        }
                    }
                }  // end for col
            }  //end for row
            Debug.WriteLine(" prepared water. Ready to route for erosion and deposition");
            all_grids = (GlobalMethods.nr) * (GlobalMethods.nc);
            memberdepressionnotconsidered = 0;
            int runner = 0;
            if (NA_anywhere_in_soil() == true) { Debug.WriteLine("NA found before sorted row col loop in water erosed"); }
            for (runner = GlobalMethods.number_of_data_cells - 1; runner >= 0; runner--)
            {     // the GlobalMethods.index is sorted from low to high values, but flow goes from high to low
                int row, col;
                if (GlobalMethods.index[runner] != -9999)
                {

                    row = GlobalMethods.row_index[runner]; col = GlobalMethods.col_index[runner];
                    //Debug.WriteLine(runner + " " + row + "  " + col + " GlobalMethods.nr " + GlobalMethods.nr + " GlobalMethods.nc " + GlobalMethods.nc + " GlobalMethods.nr*GlobalMethods.nc " + GlobalMethods.nr * GlobalMethods.nc + " data cells " + GlobalMethods.number_of_data_cells); 
                    if (GlobalMethods.t == 4 && row == 186 && col == 72) { diagnostic_mode = 1; }
                    else { diagnostic_mode = 0; }
                    powered_slope_sum = 0; max_allowed_erosion = 0; dz_min = -9999.99;
                    direct = 20; dz_max = -10; dhtemp = -99999.99; maximum_allowed_deposition = -9999.99;
                    if (GlobalMethods.depression[row, col] < 0) { GlobalMethods.depression[row, col] = 0; }
                    if ((GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 0] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 0] == col) ||
                        (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 1] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 1] == col) ||
                        (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 2] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 2] == col) ||
                        (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 3] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 3] == col) ||
                        (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 4] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 4] == col))
                    {
                        if (depressionconsidered[GlobalMethods.depression[row, col]] == 0)
                        {
                            //diagnostic_mode = 1;
                            depressionnumber = GlobalMethods.depression[row, col];
                            depressionconsidered[depressionnumber] = 1;
                            if (diagnostic_mode == 1) { Debug.WriteLine(" now considering dep " + depressionnumber + " GlobalMethods.index " + runner); }
                            update_depression(depressionnumber);
                            if (depressionsum_sediment_m == 0)
                            {
                                leave_depression_alone(depressionnumber); depressions_alone++;
                            }
                            else
                            {
                                if (depressionsum_sediment_m >= needed_to_fill_depression_m) { fill_depression(depressionnumber); depressions_filled++; }
                                else { delta_depression(depressionnumber); depressions_delta++; }
                            }
                        }
                        //all cells of this lake have now been considered, except the outlets
                    }
                    if (GlobalMethods.depression[row, col] < 0) { Debug.WriteLine(" error: negative depression value " + GlobalMethods.depression[row, col] + " at " + row + " " + col); minimaps(row, col); }
                    // this check indicates a problem with the resetting of cells involved in a delta
                    if (GlobalMethods.depression[row, col] == 0 ||
                                                (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 0] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 0] == col) ||
                                                (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 1] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 1] == col) ||
                                                (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 2] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 2] == col) ||
                                                (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 3] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 3] == col) ||
                                                (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 4] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 4] == col))
                    { //for all cells outside a depression and for outlets, we use the stream power equations based on a multiple flow (D8) template
                        //if (row == 24 && col == 81) { Debug.WriteLine(" looking around cell " + row + " " + col); minimaps(row, col); }
                        for (i = (-1); i <= 1; i++)
                        {
                            for (j = (-1); j <= 1; j++)
                            {
                                dh = 0; dhtemp = -99999.99; GlobalMethods.d_x = GlobalMethods.dx;
                                if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))  //to stay within the grid and avoid the row col cell itself
                                {
                                    // below, we calculate slope_sum for all cells either not in a depression, or being a outlet
                                    // slope_sum is needed to calculate flow in a multiple flow environment until someone thinks of something better
                                    // if (diagnostic_mode == 1) { Debug.WriteLine("checking " + (row + i) + " " + (col + j) + " from cell " + row + " " + col); }
                                    if (GlobalMethods.depression[row, col] < 0) { Debug.WriteLine(" lakes error: cell has depression < 1"); } //out_integer("wrong_lakes.asc", depression); 
                                    if (GlobalMethods.depression[row, col] == 0)
                                    {    // if the cell is not in a depression (it could be in a depression as an outlet)
                                        if (GlobalMethods.dtm[row + i, col + j] != -9999)
                                        {  //if the cell has no NODATA
                                            if (only_waterflow_checkbox.Checked)
                                            {
                                                dh = GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + i, col + j]; // in the case that we are not interested in erosion and deposition, then there is no ero and sed to query                                            }
                                            }
                                            else
                                            {
                                                dh = (GlobalMethods.dtm[row, col] + GlobalMethods.dz_ero_m[row, col] + GlobalMethods.dz_sed_m[row, col]) - (GlobalMethods.dtm[row + i, col + j] + GlobalMethods.dz_ero_m[row + i, col + j] + GlobalMethods.dz_sed_m[row + i, col + j]);    // diff @ this moment 
                                            }
                                            if (dh < 0)  // we are looking at a higher neighbour
                                            {
                                                if (dh > maximum_allowed_deposition) { maximum_allowed_deposition = dh; }   // we keep track of the minimum difference in altitude between this cell and its lowest higher neighbour - we will not raise it more, even if we would like to when the Courant criterion is violated
                                            } // end if dh
                                            if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }   // for non-cardinal neighbours, we use the adapted length
                                            if (dh > 0)
                                            {  // i j is a lower neighbour
                                                if (dh > max_allowed_erosion - dh_tol) { max_allowed_erosion = (dh - dh_tol); }  // we keep track of the minimum difference in current altitude between this cell and its highest lower neighbour - we will not erode it more, even if we would like to
                                                                                                                                 //if (diagnostic_mode == 1) { Debug.WriteLine("cell " + row + " " + col + " GlobalMethods.dtm " + GlobalMethods.dtm[row, col] + " now " + (GlobalMethods.dtm[row, col] + GlobalMethods.dz_sed_m[row, col] + GlobalMethods.dz_ero_m[row, col]) + " nb " + (row + i) + " " + (col + j) + " GlobalMethods.dtm " + GlobalMethods.dtm[row + i, col + j] + " is now " + (GlobalMethods.dtm[row + i, col + j] + GlobalMethods.dz_ero_m[row + i, col + j] + GlobalMethods.dz_sed_m[row + i, col + j])); }    
                                                dh = dh / GlobalMethods.d_x;
                                                dh = Math.Pow(dh, conv_fac);
                                                powered_slope_sum = powered_slope_sum + dh;
                                            }//end if dh  
                                        }//end if novalues
                                    }  // end if not in depression
                                    if ((GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 0] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 0] == col)
                                        || (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 1] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 1] == col)
                                        || (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 2] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 2] == col)
                                        || (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 3] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 3] == col)
                                        || (GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 4] == row && GlobalMethods.drainingoutlet_col[GlobalMethods.depression[row, col], 4] == col))
                                    {    // this cell is one of the draining outlets and is only allowed to drain to cells not in the lake																											
                                         // if the lake has been filled at this time, then all its (by now non-lake) cells have an altitude > outlet, and will not be considered for that reason
                                        if (GlobalMethods.depression[row + i, col + j] != GlobalMethods.depression[row, col])
                                        {
                                            if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                                            dh = GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + i, col + j];
                                            if (dh > 0)
                                            {// i j is a lower neighbour
                                                dh = dh / GlobalMethods.d_x;                      // dh is now equal to slope
                                                dh = Math.Pow(dh, conv_fac);            // dh is nu de helling tot de macht conv fac
                                                powered_slope_sum = powered_slope_sum + dh;
                                            } // end if lower nb
                                        } //end if nb not within depression
                                    } // end if drainingoutlet
                                }// end if boundaries
                            }//end for j
                        }//end for i, we now know slope sum for this cell. We have included cells that are in a lake in this calculation. //ArT should we replace their altitude with depressionlevel?
                        // (row == 24 && col == 81) { Debug.WriteLine("passed"); }
                        if (maximum_allowed_deposition == -9999.99) { maximum_allowed_deposition = 0; } else { maximum_allowed_deposition = -maximum_allowed_deposition; }
                        if (max_allowed_erosion < 0) { max_allowed_erosion = -dh_tol; } else { max_allowed_erosion = -max_allowed_erosion; }
                        //if (diagnostic_mode == 1) { Debug.WriteLine(" slopesum = " + slope_sum + " maximum deposition " + maximum_allowed_deposition + " maximum erosion " + max_allowed_erosion); }

                        // we are now prepared to actually calculate erosion and deposition: we can calculate how much water and sediment is redistributed using slope_sum
                        if (NA_in_soil(row, col) == true) { Debug.WriteLine("NA found before eroding " + row + " " + col); }
                        for (i = (-1); i <= 1; i++)
                        {
                            for (j = (-1); j <= 1; j++)
                            {
                                dh = 0; fraction = 0; transport_capacity_kg = 0;
                                sediment_transported = 0; detachment_rate = 0;
                                GlobalMethods.d_x = GlobalMethods.dx;
                                if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                                {  //boundaries
                                    //if (row == 24 && col == 81) { Debug.WriteLine("entered" + i + j); }
                                    if (GlobalMethods.dtm[row + i, col + j] != -9999)
                                    {
                                        if (only_waterflow_checkbox.Checked)
                                        {
                                            dh = GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + i, col + j];
                                        }
                                        else
                                        {
                                            dh = (GlobalMethods.dtm[row, col] + GlobalMethods.dz_ero_m[row, col] + GlobalMethods.dz_sed_m[row, col]) - (GlobalMethods.dtm[row + i, col + j] + GlobalMethods.dz_ero_m[row + i, col + j] + GlobalMethods.dz_sed_m[row + i, col + j]);
                                        }
                                        if (dh > 0)
                                        {  //we have found one of the lower nbs
                                            //if (row == 24 && col == 81) { Debug.WriteLine("this is a lower nb " + i + j + "dh" + dh + " " + GlobalMethods.waterflow_m3[row, col]); }
                                            if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                                            if ((GlobalMethods.depression[row, col] != 0 && GlobalMethods.depression[row + i, col + j] != GlobalMethods.depression[row, col]) || (GlobalMethods.depression[row, col] == 0))
                                            {   //if cell == outlet of current lake and nb not member of that lake OR if not a lake member


                                                // Now, we first calculate the fraction of water and sediment that goes from row, col to row+i to col+j , always using current altitudes
                                                // Then, we calculate the actual amounts of water and sediment, and with that, using the stream power equation, the transport capacity
                                                // In future, the Hjülstrom diagram can be used to give texture-dependent erosion thresholds (or selectivity)

                                                dh /= GlobalMethods.d_x;  //dh is now slope
                                                fraction = Math.Pow(dh, conv_fac) / powered_slope_sum;
                                                if (GlobalMethods.waterflow_m3[row, col] < 0) { GlobalMethods.waterflow_m3[row, col] = 0; }    // this can have happened if water enters a drier zone in the landscape
                                                flow_between_cells_m3_per_m = fraction * GlobalMethods.waterflow_m3[row, col] / GlobalMethods.dx;
                                                if (GlobalMethods.depression[row + i, col + j] == 0)
                                                {  // if receiving cell is not in a depression, its waterflow is increased 
                                                    GlobalMethods.waterflow_m3[row + i, col + j] += flow_between_cells_m3_per_m * GlobalMethods.dx;
                                                }
                                                if (GlobalMethods.depression[row + i, col + j] != 0)
                                                {  // if receiving cell is in a depression, its outlets' waterflow is increased 
                                                    currentdepression = Math.Abs(GlobalMethods.depression[row + i, col + j]); // this Abs stuff should not be necessary and is included for stability!
                                                    int outletcounter = 0;
                                                    while (GlobalMethods.drainingoutlet_col[currentdepression, outletcounter] != -1)
                                                    {
                                                        outletcounter++;
                                                        if (outletcounter == 5) { break; }
                                                    }
                                                    for (int iter = 0; iter < outletcounter; iter++) // for all outlets of this depression, divide that amount of water over them
                                                    {
                                                        GlobalMethods.waterflow_m3[GlobalMethods.drainingoutlet_col[currentdepression, iter], GlobalMethods.drainingoutlet_col[currentdepression, iter]] += GlobalMethods.dx * flow_between_cells_m3_per_m / outletcounter;
                                                    }
                                                }

                                                if (only_waterflow_checkbox.Checked == false)
                                                {

                                                    organic_in_transport = fraction * (GlobalMethods.old_SOM_in_transport_kg[row, col] + GlobalMethods.young_SOM_in_transport_kg[row, col]);    //all in kg
                                                                                                                                                                    //so far, organic in transport does not count towards the transport capacity. We can have infinite amounts of it in transport

                                                    transport_capacity_kg = advection_erodibility * (GlobalMethods.bulkdensity[row, col, 0] * GlobalMethods.dx * GlobalMethods.dx) * (Math.Pow(flow_between_cells_m3_per_m, m) * Math.Pow(dh, n)); // in a departure from literature, the erosion threshold is only evaluated if erosion actually occurs
                                                    if (transport_capacity_kg < 0)
                                                    {
                                                        transport_capacity_kg = 0;
                                                        Debug.WriteLine(" Warning: negative transport capacity at" + row + " " + col);
                                                    }  // this should never happen
                                                       // We now compare transport_capacity with the total amount of sediment in transport, to determine whether we will have erosion or deposition or nothing
                                                    total_sediment_in_transport_kg = 0;

                                                    for (size = 0; size < GlobalMethods.n_texture_classes; size++)
                                                    {
                                                        total_sediment_in_transport_kg += fraction * GlobalMethods.sediment_in_transport_kg[row, col, size];                     //all in kg
                                                    }
                                                    if (transport_capacity_kg == total_sediment_in_transport_kg)
                                                    {
                                                        // neither erosion nor deposition, simply transport
                                                        for (size = 0; size < GlobalMethods.n_texture_classes; size++)
                                                        {
                                                            GlobalMethods.sediment_in_transport_kg[row + i, col + j, size] += fraction * GlobalMethods.sediment_in_transport_kg[row, col, size];  //all in kg 
                                                        }
                                                        GlobalMethods.old_SOM_in_transport_kg[row + i, col + j] += fraction * GlobalMethods.old_SOM_in_transport_kg[row, col];  //all in kg
                                                        GlobalMethods.young_SOM_in_transport_kg[row + i, col + j] += fraction * GlobalMethods.young_SOM_in_transport_kg[row, col];  //all in kg
                                                    }
                                                    if (transport_capacity_kg > total_sediment_in_transport_kg)
                                                    {

                                                        //erosion
                                                        //in case of desired erosion, we first evaluate whether we exceed the erosion threshold
                                                        if ((transport_capacity_kg - total_sediment_in_transport_kg) > erosion_threshold_kg)
                                                        {

                                                            //first, calculate how much we are going to erode. Not as much as we want to if the soil is protected by rocks or plants
                                                            rock_fraction = GlobalMethods.texture_kg[row, col, 0, 0] / (GlobalMethods.texture_kg[row, col, 0, 0] + GlobalMethods.texture_kg[row, col, 0, 1] + GlobalMethods.texture_kg[row, col, 0, 2] + GlobalMethods.texture_kg[row, col, 0, 3] + GlobalMethods.texture_kg[row, col, 0, 4]);
                                                            if (guiVariables.Version_lux_checkbox == false)
                                                            {
                                                                mass_to_be_eroded = (transport_capacity_kg - total_sediment_in_transport_kg)
                                                                * Math.Exp(-rock_protection_constant * rock_fraction)
                                                                * Math.Exp(-bio_protection_constant * 0);
                                                            }

                                                            else
                                                            {  //for Luxemburg version, here we additially protect soil from erosion by its cover of 'bad' organic matter as litter (i.e. in top layer)

                                                                // MvdM litter fraction is determined by the total amount of litter as fraction of the mineral soil in the top layer. This might be changed, because mineral content is variable and indepent of litter quantity

                                                                //XIA change this number to 0.25 as well. For GlobalMethods.creep and no GlobalMethods.creep
                                                                double litter_characteristic_protection_mass_kg_m2 = 0.01; // based on average litter contents in Luxembourg
                                                                double litter_characteristic_protection_mass_kg = litter_characteristic_protection_mass_kg_m2 * GlobalMethods.dx * GlobalMethods.dx;
                                                                double litter_protection_fraction = Math.Exp(-litter_characteristic_protection_mass_kg / (GlobalMethods.litter_kg[row, col, 0] + GlobalMethods.litter_kg[row, col, 1]));

                                                                //double litter_fraction = (GlobalMethods.litter_kg[row, col, 0] + GlobalMethods.litter_kg[row, col, 1]) / (GlobalMethods.litter_kg[row, col, 0] + GlobalMethods.litter_kg[row, col, 1] + total_layer_mass(row, col, 0));

                                                                //double litter_fraction = (GlobalMethods.old_SOM_kg[row, col, 0] + GlobalMethods.young_SOM_kg[row, col, 0]) / total_layer_mass(row, col, 0);
                                                                //LUX Xia you have to set this parameter here in the code. Value between 0-1.
                                                                //double litter_protection_constant = 0.5;

                                                                mass_to_be_eroded = (transport_capacity_kg - total_sediment_in_transport_kg)
                                                                * Math.Exp(-rock_protection_constant * rock_fraction)
                                                                * Math.Exp(-bio_protection_constant * 0)
                                                                * Math.Exp(-litter_protection_fraction);
                                                            }
                                                            //DEV possible: - bio_protection_constant * vegetation_cover_fraction
                                                            //Debug.WriteLine("eroding " + mass_to_be_eroded + " rock exp " + Math.Exp(-rock_protection_constant * rock_fraction) + " bio exp " + Math.Exp(-bio_protection_constant * 1));
                                                            // second, calculate how the mass to be eroded is taken from the different size fractions: selectivity
                                                            // if total transport capacity is small, only the finer fractions will be eroded (selectivity with diameter to power 0.5). For larger transport capacities, selectivity decreases (diameter to power 0 = equal between fractions)
                                                            // more info in excel file in dropbox. 
                                                            double constant_b1 = 0.5 * Math.Exp(constant_selective_transcap * transport_capacity_kg);
                                                            double sum_diameter_power = 0;
                                                            for (size = 0; size < 5; size++)
                                                            {
                                                                sum_diameter_power += 1 / Math.Pow(upper_particle_size[size], constant_b1);
                                                            }
                                                            double clayeroded_0_kg = 0, claypresent_0_kg = 0, clayeroded_1_kg = 0, claypresent_1_kg = 0;
                                                            for (size = 0; size < 5; size++)
                                                            {
                                                                selectivity_fraction = (1 / Math.Pow(upper_particle_size[size], constant_b1)) / sum_diameter_power;    // unit [-]
                                                                if (GlobalMethods.texture_kg[row, col, 0, size] >= selectivity_fraction * mass_to_be_eroded)
                                                                {    // typical situation
                                                                    if (size > 2)
                                                                    {
                                                                        clayeroded_0_kg += selectivity_fraction * mass_to_be_eroded;
                                                                        claypresent_0_kg += GlobalMethods.texture_kg[row, col, 0, size];
                                                                    }
                                                                    total_mass_eroded[size] += selectivity_fraction * mass_to_be_eroded;
                                                                    GlobalMethods.texture_kg[row, col, 0, size] -= selectivity_fraction * mass_to_be_eroded;   // unit [kg]
                                                                    GlobalMethods.sediment_in_transport_kg[row + i, col + j, size] += selectivity_fraction * mass_to_be_eroded;  // unit [kg]
                                                                }
                                                                else
                                                                {    // exceptional. If we want to erode more than present in the layer, we will take it from one layer down.
                                                                     //this is to avoid exceptionally thin rocky layers blocking all erosion
                                                                     //we will then first erode everything from the top layer (layer "0") and then erode from the second layer  (i.e. layer "1").
                                                                    total_mass_eroded[size] += GlobalMethods.texture_kg[row, col, 0, size];
                                                                    double left = (selectivity_fraction * mass_to_be_eroded) - GlobalMethods.texture_kg[row, col, 0, size]; // unit [kg]
                                                                    GlobalMethods.sediment_in_transport_kg[row + i, col + j, size] += GlobalMethods.texture_kg[row, col, 0, size];
                                                                    if (size > 2)
                                                                    {
                                                                        clayeroded_0_kg += GlobalMethods.texture_kg[row, col, 0, size];
                                                                        claypresent_0_kg += 0;
                                                                    }
                                                                    GlobalMethods.texture_kg[row, col, 0, size] = 0;
                                                                    if (GlobalMethods.texture_kg[row, col, 1, size] >= left)
                                                                    {   // typical
                                                                        total_mass_eroded[size] += left;
                                                                        if (size > 2)
                                                                        {
                                                                            clayeroded_1_kg += left;
                                                                            claypresent_1_kg += GlobalMethods.texture_kg[row, col, 1, size] - left;
                                                                        }
                                                                        GlobalMethods.texture_kg[row, col, 1, size] -= left;  // unit [kg]
                                                                        GlobalMethods.sediment_in_transport_kg[row + i, col + j, size] += left;  // unit [kg]
                                                                    }
                                                                    else
                                                                    {
                                                                        total_mass_eroded[size] += GlobalMethods.texture_kg[row, col, 1, size];
                                                                        GlobalMethods.sediment_in_transport_kg[row + i, col + j, size] += GlobalMethods.texture_kg[row, col, 1, size];// unit [kg]
                                                                        if (size > 2)
                                                                        {
                                                                            clayeroded_1_kg += GlobalMethods.texture_kg[row, col, 1, size];
                                                                            claypresent_1_kg += 0;
                                                                        }
                                                                        GlobalMethods.texture_kg[row, col, 1, size] = 0;
                                                                    }
                                                                }
                                                            }

                                                            //organic matter is eroded as a fraction of total OM. That fraction equals the fraction of clay eroded from the layer
                                                            //the assumption underlying this is that clay and humus are bound in aggregates
                                                            //this does not cover: LMW SOM, peat or large woody debris
                                                            double clayerodedfraction_0 = clayeroded_0_kg / (clayeroded_0_kg + claypresent_0_kg);
                                                            double clayerodedfraction_1 = clayeroded_1_kg / (clayeroded_1_kg + claypresent_1_kg);
                                                            if (Double.IsNaN(clayerodedfraction_0))
                                                            {
                                                                clayerodedfraction_0 = 0;
                                                                //Debug.WriteLine(" this should not have happened - no OM erosion possible"); 
                                                            }
                                                            if (Double.IsNaN(clayerodedfraction_1)) { clayerodedfraction_1 = 0; }
                                                            //if (row == 62 && col == 78) { Debug.WriteLine(clayerodedfraction_0 + "  " + clayerodedfraction_1); displaysoil(row, col); }
                                                            GlobalMethods.old_SOM_in_transport_kg[row, col] += GlobalMethods.old_SOM_kg[row, col, 0] * clayerodedfraction_0 + GlobalMethods.old_SOM_kg[row, col, 1] * clayerodedfraction_1;
                                                            GlobalMethods.young_SOM_in_transport_kg[row, col] += GlobalMethods.young_SOM_kg[row, col, 0] * clayerodedfraction_0 + GlobalMethods.young_SOM_kg[row, col, 1] * clayerodedfraction_1;
                                                            total_mass_eroded[5] += GlobalMethods.old_SOM_kg[row, col, 0] * clayerodedfraction_0 + GlobalMethods.old_SOM_kg[row, col, 1] * clayerodedfraction_1;
                                                            total_mass_eroded[6] += GlobalMethods.young_SOM_kg[row, col, 0] * clayerodedfraction_0 + GlobalMethods.young_SOM_kg[row, col, 1] * clayerodedfraction_1;
                                                            GlobalMethods.old_SOM_kg[row, col, 0] *= 1 - clayerodedfraction_0;
                                                            GlobalMethods.young_SOM_kg[row, col, 0] *= 1 - clayerodedfraction_0;
                                                            GlobalMethods.old_SOM_kg[row, col, 1] *= 1 - clayerodedfraction_1;
                                                            GlobalMethods.young_SOM_kg[row, col, 1] *= 1 - clayerodedfraction_1;
                                                            //if (row == 62 && col == 78) displaysoil(row,col);
                                                        }
                                                        else
                                                        {
                                                            //do nothing. We wanted to erode, but not enough so to actually exceed the threshold and actually do that
                                                        }
                                                    }
                                                    if (transport_capacity_kg < total_sediment_in_transport_kg)
                                                    {
                                                        //deposition
                                                        //Debug.WriteLine("deposition");
                                                        //first, calculate how much we are going to keep in transport. This is the way that selectivity works now. 
                                                        double sum_diameter_power = 0, clay_deposited = 0, clay_transported = 0;
                                                        for (size = 0; size < 5; size++)
                                                        {
                                                            sum_diameter_power += 1 / Math.Pow(upper_particle_size[size], 0.5);
                                                        }
                                                        for (size = 0; size < 5; size++)
                                                        {
                                                            selectivity_fraction = (1 / Math.Pow(upper_particle_size[size], 0.5)) / sum_diameter_power;    // unit [-]
                                                            potential_transported_amount_kg = selectivity_fraction * transport_capacity_kg;                      // unit [kg]
                                                            if (potential_transported_amount_kg < GlobalMethods.sediment_in_transport_kg[row, col, size])
                                                            {
                                                                total_mass_deposited_kg[size] += GlobalMethods.sediment_in_transport_kg[row, col, size] - potential_transported_amount_kg;
                                                                GlobalMethods.texture_kg[row, col, 0, size] += GlobalMethods.sediment_in_transport_kg[row, col, size] - potential_transported_amount_kg;        // unit [kg]
                                                                GlobalMethods.sediment_in_transport_kg[row + i, col + j, 0] = potential_transported_amount_kg;                                    // unit [kg]  

                                                                if (size > 2)
                                                                {
                                                                    clay_deposited += GlobalMethods.sediment_in_transport_kg[row, col, size] - potential_transported_amount_kg;
                                                                    clay_transported += potential_transported_amount_kg;
                                                                }

                                                            }
                                                            else
                                                            {
                                                                //do nothing. We keep the sediment in transport, and do not deposit anything
                                                            }

                                                        }
                                                        // now organic matter
                                                        if (!(clay_deposited == 0 && clay_transported == 0))
                                                        {
                                                            double claydepfraction = clay_deposited / (clay_deposited + clay_transported);
                                                            total_mass_deposited_kg[5] += GlobalMethods.young_SOM_in_transport_kg[row, col] * claydepfraction;
                                                            total_mass_deposited_kg[6] += GlobalMethods.old_SOM_in_transport_kg[row, col] * claydepfraction;
                                                            GlobalMethods.young_SOM_kg[row, col, 0] += GlobalMethods.young_SOM_in_transport_kg[row, col] * claydepfraction;
                                                            GlobalMethods.old_SOM_kg[row, col, 0] += GlobalMethods.old_SOM_in_transport_kg[row, col] * claydepfraction;
                                                            GlobalMethods.young_SOM_in_transport_kg[row, col] *= 1 - claydepfraction;
                                                            GlobalMethods.old_SOM_in_transport_kg[row, col] *= 1 - claydepfraction;
                                                        }
                                                        else  //so both of the clay numbers are zero. This could be for many reasons: there is no clay to erode - there is absolutely no ero   - among them
                                                        {
                                                            //do nothing
                                                            //may be a problem if the landscape is simply clay-less: in that case, we do want to be able to erode OM.
                                                        }

                                                    }
                                                } // end if else : also erosion and deposition considered

                                            }

                                            // 4. Indien oververzadigd: depositie. Berekenen van de doorgaande massa van iedere textuurklasse, op basis van 1/d0.5 (zie Excel). 
                                            // 4b. Vergelijken van doorgaande massa met massa aanwezig in transport per textuurfractie. Indien teveel aanwezig, afwerpen. 
                                            // 4c. Organische stof afwerpen propoertioneel met de afzettingsfractie van de beide kleifracties. (Dus als er 30% van de klei in transport blijft, dan ook 30% van de OM).
                                            // Dit leidt bij de kleifractie slechts zelden tot afzetting. 

                                            // Depressies: volledige afzetting van materiaal dat in transport is. 
                                            // Instabiliteit: geen garantie dat dit niet gebeurt. Smearing kan er bij gezet worden. 
                                            // Gravelafzettingen: volgens pdf een rho van 2.7. Afgeronde gravel afzettingen van rivieren kunnen die heel laag hebben. 

                                        } //end`dH > 000
                                    }//end if novalues
                                }//end if boundaries
                            }//end for j
                        }//end for i
                        if (NA_in_soil(row, col) == true) { Debug.WriteLine("NA found after eroding " + row + " " + col); }
                        //if (row == 24 && col == 81) { Debug.WriteLine("passed"); }
                    } // end if not in a lake or a lake outlet (all other lake cells have been considered before
                } //end if nodata
            }//end for GlobalMethods.index
             // all cells have now been considered in order of (original) altitude. We must still recalculate their thicknesses and recalculate altitude. While doing that, we should count how much erosion and deposition there has been. 
            volume_eroded = 0; sediment_exported = 0; volume_deposited = 0;
            total_average_altitude = 0; total_altitude = 0;
            total_rain = 0; total_evap = 0; total_infil = 0; total_outflow = 0;
            wet_cells = 0; eroded_cells = 0; deposited_cells = 0;
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    if (GlobalMethods.dtm[row, col] != -9999)
                    {
                        if (only_waterflow_checkbox.Checked == false)
                        {
                            //erosion and deposition affect only the top two layers of soil. All others: unaffected.
                            //So, we calculate the difference between the original and final thicknesses of these two layers to calculate GlobalMethods.dz_ero_m and GlobalMethods.dz_sed_m. 
                            //We already knew how much mass was involved in ero and sed, but we need the volumes to update the GlobalMethods.dtm.

                            for (i = 0; i < 2; i++)
                            {
                                double pastlayer = GlobalMethods.layerthickness_m[row, col, i];
                                GlobalMethods.layerthickness_m[row, col, i] = thickness_calc(row, col, i);
                                if (pastlayer < GlobalMethods.layerthickness_m[row, col, i])  //if there is deposition in volume terms
                                {
                                    GlobalMethods.dz_sed_m[row, col] += GlobalMethods.layerthickness_m[row, col, i] - pastlayer;  // leading to positive values for GlobalMethods.dz_sed_m, which is what we want
                                }
                                else
                                {
                                    GlobalMethods.dz_ero_m[row, col] += GlobalMethods.layerthickness_m[row, col, i] - pastlayer;  //leading to negative values for GlobalMethods.dz_ero_m, which is what we want
                                }
                            }
                            //now GlobalMethods.dz_ero_m and GlobalMethods.dz_sed_m hold the changed altitudes. 

                            volume_eroded += GlobalMethods.dz_ero_m[row, col];
                            volume_deposited += GlobalMethods.dz_sed_m[row, col];
                            GlobalMethods.dtmchange[row, col] += GlobalMethods.dz_ero_m[row, col] + GlobalMethods.dz_sed_m[row, col];  //attention: LAKE_sed and GlobalMethods.dz_sed_m are treated differently. 
                            GlobalMethods.dtm[row, col] += GlobalMethods.dz_ero_m[row, col] + GlobalMethods.dz_sed_m[row, col];                           //No need to add lake_sed to GlobalMethods.dtm in the next line
                            GlobalMethods.sum_water_erosion[row, col] += GlobalMethods.dz_ero_m[row, col] + GlobalMethods.dz_sed_m[row, col] + GlobalMethods.lake_sed_m[row, col];

                            if (-GlobalMethods.dz_ero_m[row, col] > guiVariables.Timeseries.Timeseries_erosion_threshold) { eroded_cells++; }
                            if (GlobalMethods.dz_sed_m[row, col] + GlobalMethods.lake_sed_m[row, col] > guiVariables.Timeseries.Timeseries_deposition_threshold) { deposited_cells++; }
                        }
                        if (check_space_rain.Checked == true) { total_rain += GlobalMethods.rain[row, col]; }
                        total_rain += rain_value_m;
                        if (check_space_evap.Checked == true) { total_evap += GlobalMethods.evapotranspiration[row, col]; }
                        total_evap += evap_value_m;
                        if (check_space_infil.Checked == true) { total_infil += GlobalMethods.infil[row, col]; }
                        total_infil += infil_value_m;
                        if (GlobalMethods.waterflow_m3[row, col] * GlobalMethods.dx * GlobalMethods.dx > guiVariables.Timeseries.Timeseries_waterflow_threshold) { wet_cells++; }
                    } // end for nodata
                }   // end for col
            } // end for row
            total_rain *= GlobalMethods.dx * GlobalMethods.dx;   // m3
            total_evap *= GlobalMethods.dx * GlobalMethods.dx;   // m3
            total_infil *= GlobalMethods.dx * GlobalMethods.dx;  // m3
            total_outflow = total_rain - total_evap - total_infil;
            //Debug.WriteLine("\n--erosion and deposition overview--");
            //Debug.WriteLine("GlobalMethods.rain " + total_rain + " evap " + total_evap + " total_infil " + total_infil);
            if (only_waterflow_checkbox.Checked == false)
            {
                double total_kg_eroded = 0, total_kg_deposited = 0;
                for (size = 0; size < 5; size++)
                {
                    total_kg_eroded += total_mass_eroded[size];
                    total_kg_deposited += total_mass_deposited_kg[size];
                }

                //Debug.WriteLine(" number of dhmax erosion errors: " + + "\n" ,dhmax_errors); 
                //Debug.WriteLine(" filled " + + " of " + + " depressions, %.3f sediment used for %.3f depressionvolume\n",depressions_filled,totaldepressions,sediment_filled,depressionvolume_filled); 
                //Debug.WriteLine(" sedimented into " + + " of " + + " depressions, %.3f sediment used\n",depressions_delta,totaldepressions,sediment_delta);
                //Debug.WriteLine(" left alone " + + " of " + + " depressions",depressions_alone,totaldepressions); 
                //Debug.WriteLine(" total %6.0f cubic metres of sediment (of max %6.0f) deposited ",(sediment_deposited+sediment_delta+sediment_filled)*GlobalMethods.dx*GlobalMethods.dx,(-sediment_produced*GlobalMethods.dx*GlobalMethods.dx)); 
                /* Debug.WriteLine(" MASS BASED [kg]:");
                 Debug.WriteLine(" SDR_all " + (total_kg_eroded - total_kg_deposited) / (total_kg_eroded));
                 if (total_mass_eroded[0] != 0) { Debug.WriteLine(" SDR_coarse " + (total_mass_eroded[0] - total_mass_deposited[0]) / (total_mass_eroded[0]) + " ero " + total_mass_eroded[0] + "kg sed " + total_mass_deposited[0] + "kg"); } else { Debug.WriteLine("no coarse transport"); }
                 if (total_mass_eroded[1] != 0) { Debug.WriteLine(" SDR_sand " + (total_mass_eroded[1] - total_mass_deposited[1]) / (total_mass_eroded[1]) + " ero " + total_mass_eroded[1] + "kg sed " + total_mass_deposited[1] + "kg"); } else { Debug.WriteLine("no sand transport"); }
                 if (total_mass_eroded[2] != 0) { Debug.WriteLine(" SDR_silt " + (total_mass_eroded[2] - total_mass_deposited[2]) / (total_mass_eroded[2]) + " ero " + total_mass_eroded[2] + "kg sed " + total_mass_deposited[2] + "kg"); } else { Debug.WriteLine("no silt transport"); }
                 if (total_mass_eroded[3] != 0) { Debug.WriteLine(" SDR_clay " + (total_mass_eroded[3] - total_mass_deposited[3]) / (total_mass_eroded[3]) + " ero " + total_mass_eroded[3] + "kg sed " + total_mass_deposited[3] + "kg"); } else { Debug.WriteLine("no clay transport"); }
                 if (total_mass_eroded[4] != 0) { Debug.WriteLine(" SDR_fine_clay " + (total_mass_eroded[4] - total_mass_deposited[4]) / (total_mass_eroded[4]) + " ero " + total_mass_eroded[4] + "kg sed " + total_mass_deposited[4] + "kg"); } else { Debug.WriteLine("no fine clay transport"); }

                 Debug.WriteLine(" VOLUME BASED [m3]:");
                 Debug.WriteLine(" SDR " + (volume_eroded + volume_deposited + sediment_delta + sediment_filled) / (volume_eroded));
                 //Debug.WriteLine(" as sink : %.3f ",((-sediment_delta-sediment_filled)/sediment_produced)); 
                 //Debug.WriteLine(" as sediment : %.3f ",((-sediment_deposited)/sediment_produced)); 
                 Debug.Write(" ERO " + (volume_eroded * GlobalMethods.dx * GlobalMethods.dx) + " \n");
                 Debug.Write(" SED " + (volume_deposited * GlobalMethods.dx * GlobalMethods.dx) + " \n");
                 Debug.Write(" DEL " + (sediment_delta * GlobalMethods.dx * GlobalMethods.dx) + " \n");
                 Debug.WriteLine(" FIL " + (sediment_filled * GlobalMethods.dx * GlobalMethods.dx) + " \n");
                 */
                /*if ((volume_eroded + volume_deposited + sediment_delta + sediment_filled) / volume_eroded != 0)
                {
                    Debug.WriteLine(" ALTITUDE BASED:");
                    Debug.WriteLine(" GlobalMethods.t = " + GlobalMethods.t + " number of dhmax erosion errors: " + dhmax_errors);
                    Debug.WriteLine(" on m-basis: filled " + depressions_filled + " of " + totaldepressions + " depressions, " + sediment_filled + " sediment used for " + depressionvolume_filled + " depressionvolume");
                    Debug.WriteLine(" on m-basis: sedimented into " + depressions_delta + " of " + totaldepressions + " depressions, " + sediment_delta + "  sediment used");
                } */
            }
            this.InfoStatusPanel.Text = "calc movement has been finished";
            this.out_sed_statuspanel.Text = string.Format("sed_exp {0:F0} * 1000 m3", total_sed_export * GlobalMethods.dx * GlobalMethods.dx / 1000);


            //save timeseries_outputs
            if (guiVariables.Timeseries.Timeseries_cell_waterflow_check)
            {
                timeseries_matrix[GlobalMethods.t, timeseries_order[1]] = GlobalMethods.waterflow_m3[System.Convert.ToInt32(guiVariables.Timeseries.Timeseries_textbox_cell_row), System.Convert.ToInt32(guiVariables.Timeseries.Timeseries_textbox_cell_col)];
            }
            if (guiVariables.Timeseries.Timeseries_cell_altitude_check)
            {
                timeseries_matrix[GlobalMethods.t, timeseries_order[2]] = GlobalMethods.dtm[System.Convert.ToInt32(guiVariables.Timeseries.Timeseries_textbox_cell_row), System.Convert.ToInt32(guiVariables.Timeseries.Timeseries_textbox_cell_col)];
            }
            if (guiVariables.Timeseries.Timeseries_net_ero_check)
            {
                timeseries_matrix[GlobalMethods.t, timeseries_order[3]] = volume_eroded + volume_deposited + sediment_delta_m + sediment_filled_m;
            }
            if (guiVariables.Timeseries.Timeseries_number_dep_check)
            {
                timeseries_matrix[GlobalMethods.t, timeseries_order[4]] = deposited_cells;
            }
            if (guiVariables.Timeseries.Timeseries_number_erosion_check)
            {
                timeseries_matrix[GlobalMethods.t, timeseries_order[5]] = eroded_cells;
            }
            if (guiVariables.Timeseries.Timeseries_number_waterflow_check)
            {
                timeseries_matrix[GlobalMethods.t, timeseries_order[6]] = wet_cells;
            }
            if (guiVariables.Timeseries.Timeseries_SDR_check)
            {
                timeseries_matrix[GlobalMethods.t, timeseries_order[7]] = (volume_eroded + volume_deposited + sediment_delta_m + sediment_filled_m) / volume_eroded;
            }
            if (guiVariables.Timeseries.Timeseries_total_average_alt_check)
            {
                timeseries_matrix[GlobalMethods.t, timeseries_order[8]] = total_average_altitude;
            }
            if (guiVariables.Timeseries.Timeseries_total_dep_check)
            {
                timeseries_matrix[GlobalMethods.t, timeseries_order[9]] = volume_deposited + sediment_delta_m + sediment_filled_m;
            }
            if (guiVariables.Timeseries.Timeseries_total_ero_check)
            {
                timeseries_matrix[GlobalMethods.t, timeseries_order[10]] = -volume_eroded;
            }
            if (guiVariables.Timeseries.Timeseries_total_evap_check)
            {
                timeseries_matrix[GlobalMethods.t, timeseries_order[11]] = total_evap;
            }
            if (guiVariables.Timeseries.Timeseries_total_infil_check)
            {
                timeseries_matrix[GlobalMethods.t, timeseries_order[12]] = total_infil;
            }
            if (guiVariables.Timeseries.Timeseries_total_outflow_check)
            {
                timeseries_matrix[GlobalMethods.t, timeseries_order[13]] = total_outflow;
            }
            if (guiVariables.Timeseries.Timeseries_total_rain_check)
            {
                timeseries_matrix[GlobalMethods.t, timeseries_order[14]] = total_rain;
            }

        }

        void ini_slope()   //Initialise LS parameters   
        {
            // the soil physical / hydrological / slope stability parameters:
            //	 Transmissivity, Bulk Density,              
            //   Combined Cohesion and Internal riction.
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    //currently spatially uniform
                    GlobalMethods.T_fac[row, col] = System.Convert.ToDouble(textBox_ls_trans.Text);
                    GlobalMethods.C_fac[row, col] = System.Convert.ToDouble(textBox_ls_coh.Text);
                    GlobalMethods.bulkd[row, col] = System.Convert.ToDouble(textBox_ls_bd.Text);
                    GlobalMethods.intfr[row, col] = System.Convert.ToDouble(textBox_ls_ifr.Text);

                    // below, the old parameter values for New Zealand (spatially different) are kept
                    /*if (soilmap[row, col] == -9999)
                    { //
                        soilmap[row, col] = 0;
                    }
                    GlobalMethods.T_fac[row, col] = 15 * T_act; GlobalMethods.C_fac[row, col] = 0.2 * C_act; GlobalMethods.Cs_fac[row, col] = 10; // Defaults  15;0.2;10.0;1.8;0.7
                    GlobalMethods.bulkd[row, col] = 1.8 * bulkd_act; GlobalMethods.intfr[row, col] = 0.7 * intfr_act;
                    /*if (soilmap[row,col]==1) {  // Lone Kauri 15;0.43;12.223;1.455;0.688
                        GlobalMethods.T_fac[row,col]=a_T*T_act; GlobalMethods.C_fac[row,col]=a_coh*C_act;
                        GlobalMethods.bulkd[row,col]=a_bd*bulkd_act; GlobalMethods.intfr[row,col]=a_ifr*intfr_act;
                     }
                    if (soilmap[row,col]==2) {  // Piha       18;0.21;5.976;1.447;0.678
                        GlobalMethods.T_fac[row,col]=b_T*T_act; GlobalMethods.C_fac[row,col]=b_coh*C_act;
                        GlobalMethods.bulkd[row,col]=b_bd*bulkd_act; GlobalMethods.intfr[row,col]=b_ifr*intfr_act;
                     }
                    if (soilmap[row,col]==3) {  // Nihotupu   11;0.25;13.352;1.436;0.548
                        GlobalMethods.T_fac[row,col]=c_T*T_act; GlobalMethods.C_fac[row,col]=c_coh*C_act;
                        GlobalMethods.bulkd[row,col]=c_bd*bulkd_act; GlobalMethods.intfr[row,col]=c_ifr*intfr_act;
                     }
                     if (soilmap[row,col]==4) {  //
                        GlobalMethods.T_fac[row,col]=d_T*T_act; GlobalMethods.C_fac[row,col]=d_coh*C_act;
                        GlobalMethods.bulkd[row,col]=d_bd*bulkd_act; GlobalMethods.intfr[row,col]=d_ifr*intfr_act;
                     }
                     if (soilmap[row,col]==5) {  //
                        GlobalMethods.T_fac[row,col]=e_T*T_act; GlobalMethods.C_fac[row,col]=e_coh*C_act;
                        GlobalMethods.bulkd[row,col]=e_bd*bulkd_act; GlobalMethods.intfr[row,col]=e_ifr*intfr_act;
                     }   */
                } //for
            } //for
        }

        void calculate_critical_rain()    //Calculates Critical Steady State Rainfall for Landsliding    
        {
            // from steepest local slope, contributing area and stability parameters
            // start calculation number of contributing draining cells by multiple flow algorithm
            this.InfoStatusPanel.Text = "critical rainfall calculation";
            double beta;
            //set all start q values effective precipitation at time GlobalMethods.t
            nb_ok = 0; nb_check = 0; all_grids = 0;
            maximum_allowed_deposition = -9999; dh_tol = 0.00025;
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    GlobalMethods.camf[row, col] = 1;    // contributing area multiple flow matrix = 1
                    GlobalMethods.stslope[row, col] = 0;
                    GlobalMethods.crrain[row, col] = 0;
                }
            }

            int runner;
            for (runner = GlobalMethods.number_of_data_cells - 1; runner >= 0; runner--)
            {           // the GlobalMethods.index is sorted from low to high values, but flow goes from high to low
                int row, col;
                row = GlobalMethods.row_index[runner]; col = GlobalMethods.col_index[runner];
                // into loop for surounding grids of certain grid
                // Start first the slope_sum loop for all lower neighbour grids
                powered_slope_sum = 0; max_allowed_erosion = 0; dz_min = -9999.99;
                direct = 20; dz_max = -1; dhtemp = -99999.99; maximum_allowed_deposition = (-9999.99);

                // Repeat the loop to determine flow if all draining neighbours are known
                // but do this only once
                for (i = (-1); i <= 1; i++)
                {
                    for (j = (-1); j <= 1; j++)
                    {
                        dh = 000000; dh1 = 000; dhtemp = -99999.99; GlobalMethods.d_x = GlobalMethods.dx;
                        if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) &&   // boundaries
                             ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) &&
                       !((i == 0) && (j == 0)))
                        {
                            dh = (GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + i, col + j]);
                            if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                            if (dh < 000000)
                            {// i j is a higher neighbour
                                if (dh > dz_min) { dz_min = dh; }
                                if ((dh < 000000))
                                {// i j is a higher neighbour
                                    if (dh1 > maximum_allowed_deposition) { maximum_allowed_deposition = (dh1); }
                                }
                            }
                            if (dh > 000000)
                            {// i j is a lower neighbour
                                if ((dh > 000000))
                                {
                                    if (dh1 > max_allowed_erosion - dh_tol) { max_allowed_erosion = (dh1 - dh_tol); }
                                }
                                dh = dh / GlobalMethods.d_x;
                                if (dh > dz_max) { dz_max = dh; direct = (i * 3 + 5 + j); }
                                dh = Math.Pow(dh, conv_fac);
                                powered_slope_sum = powered_slope_sum + dh;
                            }//end if
                        }//end if
                    }//end for
                }//end for
                if (maximum_allowed_deposition == -9999.99) { maximum_allowed_deposition = 0; } else { maximum_allowed_deposition = -maximum_allowed_deposition; }
                if (max_allowed_erosion == 0.0) { max_allowed_erosion = -dh_tol; } else { max_allowed_erosion = -max_allowed_erosion; }
                for (i = (-1); i <= 1; i++)
                {
                    for (j = (-1); j <= 1; j++)
                    {
                        dh = 000000; fraction = 0;
                        frac_dis = 0;
                        GlobalMethods.d_x = GlobalMethods.dx;
                        if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                        {
                            dh = (GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + i, col + j]);
                            // Multiple Flow: If there are lower neighbours start evaluating
                            if (dh > 000000)
                            { // multiple flow
                              // fraction of discharge into a neighbour grid
                                if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                                Slope = dh / GlobalMethods.d_x;
                                dh = dh / GlobalMethods.d_x;
                                dh = Math.Pow(dh, conv_fac);
                                fraction = (dh / powered_slope_sum); // multiple fow
                                frac_dis = (GlobalMethods.camf[row, col] * fraction);
                                GlobalMethods.camf[row + i, col + j] += frac_dis;
                            }//end if
                        }//end if boarders
                    }//end for j
                }//end for i
            }   // end for

            // Calculation of steepest descent local slope, 8 cell window
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    direct = 20; dz_max = -1;
                    for (i = (-1); i <= 1; i++)
                    {
                        for (j = (-1); j <= 1; j++)
                        {
                            dh = 000000;
                            if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                            {
                                dh = (GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + i, col + j]);
                                if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                                if (dh > 000000)
                                {// i j is a lower neighbour
                                    dh = dh / GlobalMethods.d_x;
                                    if (dh > dz_max) { dz_max = dh; direct = (i * 3 + 5 + j); }
                                }//end if
                            }//end if
                        }//end for
                    }//end for
                    for (i = (-1); i <= 1; i++)
                    {
                        for (j = (-1); j <= 1; j++)
                        {
                            dh = 000000;
                            if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                            {
                                dh = (GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + i, col + j]);
                                if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                                if ((i * 3 + 5 + j) == direct)
                                { // steepest descent
                                    GlobalMethods.stslope[row, col] = Math.Atan(dh / GlobalMethods.d_x);
                                    // Calculation of CRITICAL RAINFALL value = relative landslide hazard, along steepest descent local slope
                                    beta = (GlobalMethods.T_fac[row, col] * (Math.Sin(GlobalMethods.stslope[row, col])) * (GlobalMethods.dx / (GlobalMethods.camf[row, col] * GlobalMethods.dx * GlobalMethods.dx)) * GlobalMethods.bulkd[row, col] * (1 - ((Math.Sin(GlobalMethods.stslope[row, col]) - GlobalMethods.C_fac[row, col]) / ((Math.Tan(GlobalMethods.intfr[row, col]) * Math.Cos(GlobalMethods.stslope[row, col])))))); // 'valid' critical rainfall value
                                    if (Math.Tan(GlobalMethods.stslope[row, col]) > (Math.Tan(GlobalMethods.intfr[row, col]) + (GlobalMethods.C_fac[row, col] / Math.Cos(GlobalMethods.stslope[row, col])))) { beta = -99; } //unconditionally unstable
                                    if (((GlobalMethods.bulkd[row, col] * Math.Sin(GlobalMethods.stslope[row, col])) + ((1 - GlobalMethods.bulkd[row, col]) * Math.Cos(GlobalMethods.stslope[row, col]) * Math.Tan(GlobalMethods.intfr[row, col]))) <= ((GlobalMethods.bulkd[row, col]) * (GlobalMethods.C_fac[row, col]))) { beta = 99; } // unconditionally stable
                                    GlobalMethods.crrain[row, col] = (beta);
                                    //Debug.WriteLine( "critical GlobalMethods.rain for " + row + " " + col + " " + GlobalMethods.crrain[row,col] + " GlobalMethods.T_fac " + GlobalMethods.T_fac[row, col] + " GlobalMethods.stslope_sin " + Math.Sin(GlobalMethods.stslope[row, col]) + " upstream " + GlobalMethods.camf[row,col] + "\n GlobalMethods.bulkd " + GlobalMethods.bulkd[row, col] + " GlobalMethods.C_fac " + GlobalMethods.C_fac[row, col] + " GlobalMethods.intfr " + Math.Tan(GlobalMethods.intfr[row, col]) + " GlobalMethods.stslope_cos " + Math.Cos(GlobalMethods.stslope[row, col]) );
                                }
                            }//end if
                        }//end for
                    }//end for
                } // end for
            } // end for 
            GlobalMethods.out_double("critrain.asc", GlobalMethods.crrain);
        }

        void steepdesc(int rowst, int colst)
        {
            int trow;
            int tcol;
            trow = rowst;
            tcol = colst;
            xrow = 0; xcol = 0;
            powered_slope_sum = 0;
            for (i = (-1); i <= 1; i++)
            {
                for (j = (-1); j <= 1; j++)
                {
                    dh = 000000; dh1 = 000; dhtemp = -99999.99; GlobalMethods.d_x = GlobalMethods.dx;
                    if (((trow + i) >= 0) && ((trow + i) < GlobalMethods.nr) &&   // boundaries
                        ((tcol + j) >= 0) && ((tcol + j) < GlobalMethods.nc) &&
                        !((i == 0) && (j == 0)))
                    {
                        dh = (GlobalMethods.dtm[trow, tcol] - GlobalMethods.dtm[trow + i, tcol + j]);
                        if ((trow != trow + i) && (tcol != tcol + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                        if (dh < 000000)
                        {// i j is a higher neighbour
                            if (dh > dz_min) { dz_min = dh; }
                            if ((dh < 000000))
                            {// i j is a higher neighbour
                                if (dh1 > maximum_allowed_deposition) { maximum_allowed_deposition = (dh1); }
                            }
                        }
                        if (dh > 000000)
                        {// i j is a lower neighbour
                            if ((dh > 000000))
                            {
                                if (dh1 > max_allowed_erosion - dh_tol) { max_allowed_erosion = (dh1 - dh_tol); }
                            }
                            dh = dh / GlobalMethods.d_x;
                            if (dh > dz_max) { dz_max = dh; direct = (i * 3 + 5 + j); }
                            dh = Math.Pow(dh, conv_fac);
                            powered_slope_sum = powered_slope_sum + dh;
                        }//end if
                    }//end if
                }//end for
            }//end for
            if (maximum_allowed_deposition == -9999.99) { maximum_allowed_deposition = 0; } else { maximum_allowed_deposition = (maximum_allowed_deposition * (-1)); }
            if (max_allowed_erosion == 0) { max_allowed_erosion = dh_tol * -1; } else { max_allowed_erosion = (max_allowed_erosion * (-1)); }
            for (i = (-1); i <= 1; i++)
            {
                for (j = (-1); j <= 1; j++)
                {
                    dh = 000000; fraction = 0;
                    frac_dis = 0;
                    if (((trow + i) >= 0) && ((trow + i) < GlobalMethods.nr) &&   // boundaries
                           ((tcol + j) >= 0) && ((tcol + j) < GlobalMethods.nc) &&
                          !((i == 0) && (j == 0)))
                    {
                        dh = (GlobalMethods.dtm[trow, tcol] - GlobalMethods.dtm[trow + i, tcol + j]);
                        // Steepest descent only one neighbour
                        if ((i * 3 + 5 + j) == direct)
                        { //steepest descent
                            xrow = trow + i;
                            xcol = tcol + j;
                        }//end if
                    }//end if borders
                }//end for j
            }//end for i
        }

        void calculate_slide()
        {
            try
            {
                this.InfoStatusPanel.Text = "landslide calculation";
                int tell;
                //set all start q values effective precipitation at time GlobalMethods.t
                nb_ok = 0; nb_check = 0; all_grids = 0.0;
                maximum_allowed_deposition = -9999.0; dh_tol = 0.00025; erotot = 0.0;
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        GlobalMethods.slidemap[row, col] -= 1;  // terug opbouwen van 'landslide potential' bij meerdere tijdstappen
                        if (GlobalMethods.slidemap[row, col] < 0) { GlobalMethods.slidemap[row, col] = 0; }
                        GlobalMethods.ero_slid[row, col] = 0.0;
                        GlobalMethods.sed_slid[row, col] = 0.0;
                        GlobalMethods.cel_dist[row, col] = 0.0;
                        GlobalMethods.dh_slid[row, col] = 0.0;
                        GlobalMethods.sed_bud[row, col] = 0.0;
                    }
                }

                // into while loop for all grids if not all neighbours are processed
                int runner;
                for (runner = GlobalMethods.number_of_data_cells - 1; runner >= 0; runner--)
                {           // the GlobalMethods.index is sorted from low to high values, but flow goes from high to low
                    int row, col;
                    row = GlobalMethods.row_index[runner]; col = GlobalMethods.col_index[runner];

                    // into loop for surrounding grids of certain grid
                    // Start first the slope_sum loop for all lower neighbour grids
                    powered_slope_sum = 0.0; max_allowed_erosion = 0.0; dz_min = -9999.99; GlobalMethods.d_x = GlobalMethods.dx;
                    direct = 20; dz_max = -1.0; dhtemp = -99999.99; maximum_allowed_deposition = (-9999.99);
                    // Repeat the loop to determine if all neigbours are processed
                    nb_ok = 1;
                    for (i = (-1); i <= 1; i++)
                    {
                        for (j = (-1); j <= 1; j++)
                        {
                            dh = 0.000;
                            if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                            {
                                dh = (GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + i, col + j]);
                            }//end if
                        }//end for
                    }//end for
                     // Repeat the loop to determine flow if all draining neighbours are known
                     // but do this only once
                     // First loop to process slide erosion with a slope limit and steepest descent


                    slide_tot = 0.0;
                    dh_tot = 0.0;
                    steepdesc(row, col);
                    dh = (GlobalMethods.dtm[row, col] - GlobalMethods.dtm[xrow, xcol]);
                    if ((row != xrow) && (col != xcol)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                    Slope = dh / GlobalMethods.d_x;
                    if ((Slope > 0.176327) && (GlobalMethods.slidemap[row, col] < 1))
                    { // FACTOR 1 and not slided yet
                        if (GlobalMethods.watsh[row, col] == 1)
                        {
                            if ((GlobalMethods.crrain[row, col] > 0.0) && (GlobalMethods.crrain[row, col] < 0.02))
                            { // FACTOR 4 RELATIVE RISK FOR 'GOING' SET AT 0.02 m/d !!! = SCENARIO
                                if (GlobalMethods.ero_slid[row, col] > -((GlobalMethods.bulkd[row, col] * 9.81 * Math.Cos(Slope) * (Math.Tan(Slope) - Math.Tan(0.176327))) / GlobalMethods.Cs_fac[row, col]))
                                { // FACTOR 1 maximal erosion applied if more than one slide
                                    GlobalMethods.ero_slid[row, col] = -((GlobalMethods.bulkd[row, col] * 9.81 * Math.Cos(Slope) * (Math.Tan(Slope) - Math.Tan(0.176327))) / GlobalMethods.Cs_fac[row, col]); // FACTOR 1
                                    slide_tot += -((GlobalMethods.bulkd[row, col] * 9.81 * Math.Cos(Slope) * (Math.Tan(Slope) - Math.Tan(0.176327))) / GlobalMethods.Cs_fac[row, col]); // FACTOR 1
                                    dh_tot += dh;
                                    //getch();
                                }
                                while (Slope > 0.176327)
                                {   // FACTOR 1
                                    xxrow = xrow; xxcol = xcol;
                                    steepdesc(xrow, xcol);
                                    dh = (GlobalMethods.dtm[xxrow, xxcol] - GlobalMethods.dtm[xrow, xcol]);
                                    if ((xxrow != xrow) && (xxcol != xcol)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                                    Slope = dh / GlobalMethods.d_x;
                                    if (Slope > 0.176327)
                                    {// FACTOR 1 slide keeps eroding if > 10 degrees steepest descent is encountered
                                        if ((GlobalMethods.ero_slid[xxrow, xxcol] == 0.0) && (GlobalMethods.slidemap[row, col] < 1))
                                        { // has not been processed (eroded) earlier
                                            GlobalMethods.ero_slid[xxrow, xxcol] = -((GlobalMethods.bulkd[row, col] * 9.81 * Math.Cos(Slope) * (Math.Tan(Slope) - Math.Tan(0.176327))) / GlobalMethods.Cs_fac[row, col]); //FACTOR 1
                                            slide_tot += -((GlobalMethods.bulkd[row, col] * 9.81 * Math.Cos(Slope) * (Math.Tan(Slope) - Math.Tan(0.176327))) / GlobalMethods.Cs_fac[row, col]); // FACTOR 1
                                            dh_tot += dh;
                                        }
                                    }
                                    else { Slope = 0.0; GlobalMethods.sed_bud[xxrow, xxcol] += (slide_tot * -1.0); GlobalMethods.dh_slid[xxrow, xxcol] += (dh_tot); erotot += slide_tot; }
                                }//end while
                            }//end if
                        }
                    }

                }       // end for all sorted cells
                        //2 Second while loop to process slide deposition with a 'cell distance' and 'double' multiple flow
                nb_ok = 0; nb_check = 0; all_grids = 0.0; tell = 0;
                maximum_allowed_deposition = -9999.0; dh_tol = 0.00025; sedtot = 0.0; strsed = 0.0; startsed = 0.0;
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        GlobalMethods.cel_dist[row, col] = ((0.4 * GlobalMethods.dh_slid[row, col]) / GlobalMethods.dx); // FACTOR 2 calculate 'celdistance', empirical fraction of runout set at 0.4 (Lit.)
                        startsed += GlobalMethods.sed_bud[row, col]; // 'startsed'-counter = only to display initial sediment budget in ero-sed balance in model run
                    }
                }
                //2 into while loop for all grids if not all neighbours are processed
                for (runner = GlobalMethods.number_of_data_cells - 1; runner >= 0; runner--)
                {           // the GlobalMethods.index is sorted from low to high values, but flow goes from high to low
                    int row, col;
                    row = GlobalMethods.row_index[runner]; col = GlobalMethods.col_index[runner];
                    //2 into loop for surounding grids of certain grid
                    //2 Start first the slope_sum loop for all lower neighbour grids
                    powered_slope_sum = 0.0; max_allowed_erosion = 0.0; dz_min = -9999.99; GlobalMethods.d_x = GlobalMethods.dx;
                    direct = 20; dz_max = -1.0; dhtemp = -99999.99; maximum_allowed_deposition = (-9999.99);
                    nb_ok = 1;
                    for (i = (-1); i <= 1; i++)
                    {
                        for (j = (-1); j <= 1; j++)
                        {
                            dh = 0.000;
                            if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                            {
                                dh = (GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + i, col + j]);
                            }//end if
                        }//end for
                    }//end for
                     //2 Repeat the loop to determine flow if all draining neighbours are known
                     //2 but do this only once
                    if ((GlobalMethods.sed_bud[row, col] > 0.0) && (GlobalMethods.cel_dist[row, col] > 0.0))
                    {
                        if (GlobalMethods.sed_bud[row, col] < 0.00001) tell++;
                        for (i = (-1); i <= 1; i++)
                        {
                            for (j = (-1); j <= 1; j++)
                            {
                                dh = 0.000000; dh1 = 0.000; dhtemp = -99999.99; GlobalMethods.d_x = GlobalMethods.dx;
                                if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) &&   // boundaries
                                     ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) &&
                            !((i == 0) && (j == 0)))
                                {
                                    dh = (GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + i, col + j]);
                                    if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                                    if (dh < 0.000000)
                                    {// i j is a higher neighbour
                                        if (dh > dz_min) { dz_min = dh; }
                                        if ((dh < 0.000000))
                                        {// i j is a higher neighbour
                                            if (dh1 > maximum_allowed_deposition) { maximum_allowed_deposition = (dh1); }
                                        }
                                    }
                                    if (dh > 0.000000)
                                    {// i j is a lower neighbour
                                        if ((dh > 0.000000))
                                        {
                                            if (dh1 > max_allowed_erosion - dh_tol) { max_allowed_erosion = (dh1 - dh_tol); }
                                        }
                                        dh = dh / GlobalMethods.d_x;
                                        if (dh > dz_max) { dz_max = dh; direct = (i * 3 + 5 + j); }
                                        dh = Math.Pow(dh, conv_fac);
                                        powered_slope_sum = powered_slope_sum + dh;
                                    }//end if
                                }//end if
                            }//end for
                        }//end for
                        if (maximum_allowed_deposition == -9999.99) { maximum_allowed_deposition = 0.0; } else { maximum_allowed_deposition = (maximum_allowed_deposition * (-1.0)); }
                        if (max_allowed_erosion == 0.0) { max_allowed_erosion = dh_tol * -1.0; } else { max_allowed_erosion = (max_allowed_erosion * (-1.0)); }
                        for (i = (-1); i <= 1; i++)
                        {
                            for (j = (-1); j <= 1; j++)
                            {
                                dh = 0.000000; fraction = 0.0;
                                frac_dis = 0.0; frac_bud = 0.0;
                                GlobalMethods.d_x = GlobalMethods.dx;
                                if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                                {
                                    dh = (GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + i, col + j]);
                                    //Multiple Flow: If there are lower neighbours start evaluating
                                    if (dh > 0.000000)
                                    {// && (GlobalMethods.cel_dist[row,col]>0.0)) { // multiple flow, 'celdistance'
                                     // fraction of discharge into a neighbour grid
                                        if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                                        Slope = dh / GlobalMethods.d_x;
                                        dh = dh / GlobalMethods.d_x;
                                        dh = Math.Pow(dh, conv_fac);
                                        fraction = (dh / powered_slope_sum); // multiple fow

                                        if (GlobalMethods.cel_dist[row, col] <= 1.0)
                                        {
                                            frac_bud = (GlobalMethods.sed_bud[row, col] * fraction);
                                        }
                                        else
                                        {
                                            frac_bud = ((GlobalMethods.sed_bud[row, col] / GlobalMethods.cel_dist[row, col]) * fraction);
                                        }
                                        GlobalMethods.sed_bud[row + i, col + j] += ((GlobalMethods.sed_bud[row, col] * fraction) - frac_bud);
                                        GlobalMethods.sed_slid[row, col] += frac_bud;
                                        sedtot += frac_bud;
                                        if ((GlobalMethods.cel_dist[row, col] - 1.0) > 0.0)
                                        {
                                            if ((GlobalMethods.sed_bud[row + i, col + j] > 0.0) && (GlobalMethods.cel_dist[row + i, col + j] < (GlobalMethods.cel_dist[row, col] - 1.0)))
                                            {
                                                GlobalMethods.cel_dist[row + i, col + j] = (GlobalMethods.cel_dist[row, col] - 1.0);
                                            }
                                            else { GlobalMethods.cel_dist[row + i, col + j] += 0.0; }
                                        }
                                        else { GlobalMethods.cel_dist[row + i, col + j] += 0.0; }
                                        if ((GlobalMethods.camf[row, col] >= 500.0) && (GlobalMethods.sed_slid[row, col] > 0.0))
                                        { // FACTOR 3
                                            strsed += GlobalMethods.sed_slid[row, col];
                                        }
                                    }//end if
                                }//end borders
                            }//end for j
                        }//end for i
                    }//2 end if GlobalMethods.sed_bud
                } //2 end for all cells 2

                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (GlobalMethods.ero_slid[row, col] < 0.0)
                        {
                            GlobalMethods.slidemap[row, col] = 4;//potentieel terug opbouwen in # tijdstappen..., sediment kan eventueel al direct terug 'aan' gezet worden, zie hierboven **
                        }
                    }
                }
                Debug.WriteLine("Balance ero: %8.4f sed: %8.4f start:%8.4f strsed:%8.4f", erotot, sedtot, startsed, strsed);
            }
            catch
            {
                Debug.WriteLine("err_sli1");
            }

        } // end calc_slide()      

        private void calculate_tillage()
        {
            try
            {
                double mass_before = total_catchment_mass();
                this.InfoStatusPanel.Text = "tillage calculation";
                int row, col, i, j;
                double slope_sum, dz_min, d_x, dz_max, dh, fraction, temptill, tempdep, temptill_kg,
                            slope;

                nb_ok = 0; nb_check = 0;
                for (row = 0; row < GlobalMethods.nr; row++)
                {
                    for (col = 0; col < GlobalMethods.nc; col++)
                    {
                        GlobalMethods.till_result[row, col] = 0;
                        // if (GlobalMethods.dtm[row, col] < -9900 && GlobalMethods.dtm[row, col] != -9999) { Debug.WriteLine(" Cell " + row + " " + col + " has altitude " + GlobalMethods.dtm[row, col] + " till " + GlobalMethods.till_result[row, col]); }
                    }
                }

                int runner = 0;
                for (runner = GlobalMethods.number_of_data_cells - 1; runner >= 0; runner--)
                {           // the GlobalMethods.index is sorted from low to high values, but flow goes from high to low
                    row = GlobalMethods.row_index[runner]; col = GlobalMethods.col_index[runner];
                    // Debug.WriteLine("till1");

                    if (GlobalMethods.tillfields[row, col] == 1)
                    {
                        if (check_negative_weight(row, col) == true) { MessageBox.Show("negative weight in GlobalMethods.t " + GlobalMethods.t + ", row " + row + ", col " + col + ", step 1"); }

                        // 1. Mixing of the topsoil. 
                        double mixeddepth = 0, completelayerdepth = 0, newdepth = 0;
                        int completelayers = -1;

                        while (mixeddepth <= plough_depth)
                        {
                            completelayers++;
                            mixeddepth += GlobalMethods.layerthickness_m[row, col, completelayers];
                            // OSL_age[row, col, completelayers] = 0;


                        }// this will lead to incorporation of the (partial) layer below tillage horizon in completelayers parameter. So the highest number indicates the partial layer 
                         // Debug.WriteLine("till2");
                        double[] tilled_text = new double[5]; // includes soil 
                        double[] tilled_om = new double[2]; // includes OM
                        double[] alldepths = new double[completelayers]; // contains thicknesses of all layers
                        double[] fraction_mixed = new double[completelayers + 1];

                        // add material from complete layers
                        double mass_soil_before = total_soil_mass(row, col);
                        for (int lay = 0; lay < completelayers; lay++) // accounted for partial layer, only select complete layers
                        {
                            completelayerdepth += GlobalMethods.layerthickness_m[row, col, lay];
                            alldepths[lay] = GlobalMethods.layerthickness_m[row, col, lay];
                            fraction_mixed[lay] = 1;

                            for (int tex = 0; tex < 5; tex++)
                            {
                                tilled_text[tex] += GlobalMethods.texture_kg[row, col, lay, tex];
                            }
                            tilled_om[0] += GlobalMethods.old_SOM_kg[row, col, lay];
                            tilled_om[1] += GlobalMethods.young_SOM_kg[row, col, lay];

                        }
                        // Debug.WriteLine("till3");
                        // add material from partial layer and appoint mixed material, and give back material at the same time
                        double frac_ap = (plough_depth - completelayerdepth) / GlobalMethods.layerthickness_m[row, col, completelayers];
                        fraction_mixed[completelayers] = frac_ap;
                        if (frac_ap > 1)
                        {
                            Debug.WriteLine("err_ti1");
                        }
                        for (int tex = 0; tex < 5; tex++) // add partial mass of partial layer
                        {
                            tilled_text[tex] += GlobalMethods.texture_kg[row, col, completelayers, tex] * frac_ap; // add fraction from partial layer
                            GlobalMethods.texture_kg[row, col, completelayers, tex] *= (1 - frac_ap); // subtract mixed part
                            GlobalMethods.texture_kg[row, col, completelayers, tex] += tilled_text[tex] * (GlobalMethods.layerthickness_m[row, col, completelayers] * frac_ap) / plough_depth; // add part from mixed 
                        }

                        tilled_om[0] += GlobalMethods.old_SOM_kg[row, col, completelayers] * frac_ap;
                        GlobalMethods.old_SOM_kg[row, col, completelayers] *= (1 - frac_ap);
                        GlobalMethods.old_SOM_kg[row, col, completelayers] += tilled_om[0] * (GlobalMethods.layerthickness_m[row, col, completelayers] * frac_ap) / plough_depth;

                        tilled_om[1] += GlobalMethods.young_SOM_kg[row, col, completelayers] * frac_ap;
                        GlobalMethods.young_SOM_kg[row, col, completelayers] *= (1 - frac_ap);
                        GlobalMethods.young_SOM_kg[row, col, completelayers] += tilled_om[1] * (GlobalMethods.layerthickness_m[row, col, completelayers] * frac_ap) / plough_depth;

                        // Debug.WriteLine("till4");
                        for (int lay = 0; lay < completelayers; lay++)
                        {
                            for (int tex = 0; tex < 5; tex++)
                            {
                                GlobalMethods.texture_kg[row, col, lay, tex] = tilled_text[tex] * (alldepths[lay] / plough_depth);

                            }
                            GlobalMethods.old_SOM_kg[row, col, lay] = tilled_om[0] * (alldepths[lay] / plough_depth);
                            GlobalMethods.young_SOM_kg[row, col, lay] = tilled_om[1] * (alldepths[lay] / plough_depth);


                            GlobalMethods.layerthickness_m[row, col, lay] = thickness_calc(row, col, lay);
                            GlobalMethods.layerthickness_m[row, col, lay] = thickness_calc(row, col, lay);
                            newdepth += GlobalMethods.layerthickness_m[row, col, lay];
                        }
                        newdepth += GlobalMethods.layerthickness_m[row, col, completelayers];





                        // Debug.WriteLine("till5");
                        // 2. Calculate redistribution of material
                        // 2.a First calculate slope_sum for multiple flow, and remember how much lower the !currently! lowest lower neighbour is
                        slope_sum = 0; GlobalMethods.d_x = GlobalMethods.dx; dhtemp = -99999.99; nb_ok = 1; dz_max = 0; dz_min = -9999;
                        for (i = (-1); i <= 1; i++)
                        {
                            for (j = (-1); j <= 1; j++)
                            {
                                dh = 0; dhtemp = -99999.99; GlobalMethods.d_x = GlobalMethods.dx;
                                if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                                {    // boundaries
                                    if (GlobalMethods.dtm[row + i, col + j] != -9999)
                                    {
                                        dh = (GlobalMethods.dtm[row, col] + GlobalMethods.till_result[row, col] - GlobalMethods.dtm[row + i, col + j] + GlobalMethods.till_result[row + i, col + j]);
                                        if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                                        if (dh > 0)
                                        {           // i j is a lower neighbour
                                            if (dh > dz_max) { dz_max = dh; }
                                            dh = dh / GlobalMethods.d_x;
                                            dh = Math.Pow(dh, conv_fac);
                                            slope_sum = slope_sum + dh;
                                        }//end if
                                    }//end if novalues
                                }// end if boundaries
                            }//end for
                        }//end for

                        // 2.b knowing slope_sum, we can now calculate which fraction of the tilled amount goes where, and how much that is. 
                        // knowing the lowest lower neighbour of row,col lets us limit the tillage-erosion to avoid row,col becoming lower 
                        // than its lowest lower neighbour (avoiding sinks).
                        // we are also going to limit the tilled amount to avoid row+i, col+j becoming higher than its own lowest higher nb.
                        // that avoids sinks as well.
                        double mass_soil_after = total_soil_mass(row, col);
                        if (Math.Abs(mass_soil_before - mass_soil_after) > 0.0001)
                        {
                            Debug.WriteLine("err_ti2");
                        }
                        // Debug.WriteLine("till6");
                        for (i = (-1); i <= 1; i++)
                        {
                            for (j = (-1); j <= 1; j++)
                            {
                                dh = 0; fraction = 0.0;
                                frac_dis = 0.0;
                                GlobalMethods.d_x = GlobalMethods.dx;
                                if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0))) // boundaries
                                {
                                    if ((GlobalMethods.dtm)[row + i, col + j] != (-9999))
                                    {
                                        dh = (GlobalMethods.dtm[row, col] + GlobalMethods.till_result[row, col] - GlobalMethods.dtm[row + i, col + j] + GlobalMethods.till_result[row + i, col + j]);
                                        if (dh > 0.000000) // i j is a lower neighbour to which we would like to till a certain amount.
                                        {
                                            // Calculate fraction of discharge into this cell
                                            if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                                            slope = dh / GlobalMethods.d_x;
                                            dh = dh / GlobalMethods.d_x;
                                            dh = Math.Pow(dh, conv_fac);
                                            fraction = (dh / slope_sum);
                                            // Tillage erosion calculation
                                            temptill = fraction * (tilc * slope * plough_depth) * dt;    // temptill is what we would like to till from r,c to r+i,c+j
                                                                                                         // Tillage erosion correction through calculating maximum tillage: tempdep
                                            tempdep = GlobalMethods.soildepth_m[row, col];
                                            //if there is more soil than the difference between the donor cell and its lowest lower nb, limit tillage to that difference.
                                            if (tempdep > dz_max) { tempdep = dz_max; }
                                            //if there is not enough space in the receiver cell because its currently lowest higher neighbour is not high enough, 
                                            //then limit tillage to that amount. First, calculate the current altitude difference with the lowest higher neighbour of r+i,c+j 
                                            dz_min = 9999;
                                            for (alpha = (-1); alpha <= 1; alpha++)
                                            {
                                                for (beta = (-1); beta <= 1; beta++)
                                                {
                                                    if (((row + i + alpha) >= 0) && ((row + i + alpha) < GlobalMethods.nr) && ((col + j + beta) >= 0) && ((col + j + beta) < GlobalMethods.nc) && !((alpha == 0) && (beta == 0))) // boundaries
                                                    {
                                                        if (GlobalMethods.dtm[row + i + alpha, col + j + beta] != -9999)
                                                        {
                                                            if (GlobalMethods.dtm[row + i + alpha, col + j + beta] + GlobalMethods.till_result[row + i + alpha, col + j + beta] > (GlobalMethods.dtm[row + i, col + j] + GlobalMethods.till_result[row + i, col + j]))
                                                            { // we are looking at a higher neighbour of the receiver cell
                                                                if (GlobalMethods.dtm[row + i + alpha, col + j + beta] + GlobalMethods.till_result[row + i + alpha, col + j + beta] - GlobalMethods.dtm[row + i, col + j] + GlobalMethods.till_result[row + i, col + j] < dz_min)
                                                                {
                                                                    dz_min = GlobalMethods.dtm[row + i + alpha, col + j + beta] + GlobalMethods.till_result[row + i + alpha, col + j + beta] - GlobalMethods.dtm[row + i, col + j] + GlobalMethods.till_result[row + i, col + j];
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            // knowing the maximum tillage that the receiver cell can receive without blocking its own higher nbs, we limit the maximum tillage to that amount
                                            if (tempdep > dz_min) { tempdep = dz_min; }
                                            if (dz_min == 9999) { tempdep = 0.0; }     // if the receiver cell does not have higher nbs, we cannot till at all.
                                            if (tempdep < 0.0) tempdep = 0.0;
                                            if (tempdep > epsilon) { tempdep -= epsilon; }  // finally, always till just a bit less than the max allowed to prevent flat areas
                                            if (temptill > tempdep) { temptill = tempdep; } // if we want to till more than the maximum possible, we only till the maximum possible.
                                                                                            // update the corresponding grids 

                                            GlobalMethods.till_result[row, col] -= temptill;
                                            GlobalMethods.till_result[row + i, col + j] += temptill;
                                            //GlobalMethods.soildepth_m[row, col] -= temptill;
                                            //GlobalMethods.soildepth_m[row + i, col + j] += temptill;
                                            //if (GlobalMethods.soildepth_m[row, col] < 0) { GlobalMethods.soildepth_m[row, col] = 0; }
                                            //if (GlobalMethods.soildepth_m[row + i, col + j] < 0) { GlobalMethods.soildepth_m[row + i, col + j] = 0; }

                                            if (check_negative_weight(row, col) == true) { MessageBox.Show("negative weight in GlobalMethods.t " + GlobalMethods.t + ", row " + row + ", col " + col + ", step 2"); }


                                            //double dz_till_m = temptill;
                                            // Debug.WriteLine("till7");
                                            // 2.c update soil properties which are tilled
                                            // top layers are mixed, so it doesn'GlobalMethods.t matter where eroded material comes from.
                                            // problems can arise when eroded depth is larger than plough depth. 
                                            // development needed for layers with varying bulk density, in the case this occurs in an Ap horizon
                                            double mass_partial_layer, frac_eroded, total_mass_start, total_mass_end;

                                            //total_mass_start = total_soil_mass(row, col);
                                            int layero = 0;
                                            double temptill0 = temptill;
                                            while (temptill >= GlobalMethods.layerthickness_m[row, col, layero] | layero >= GlobalMethods.max_soil_layers) // hele laag wordt verwijderd, al het materiaal naar de volgende cel
                                            {
                                                for (int tex = 0; tex < 5; tex++)
                                                {
                                                    GlobalMethods.texture_kg[row + i, col + j, 0, tex] += GlobalMethods.texture_kg[row, col, layero, tex];
                                                    GlobalMethods.texture_kg[row, col, layero, tex] = 0;
                                                }
                                                GlobalMethods.young_SOM_kg[row + i, col + j, 0] += GlobalMethods.young_SOM_kg[row, col, layero];
                                                GlobalMethods.young_SOM_kg[row, col, layero] = 0;
                                                GlobalMethods.old_SOM_kg[row + i, col + j, 0] += GlobalMethods.old_SOM_kg[row, col, layero];
                                                GlobalMethods.old_SOM_kg[row, col, layero] = 0;

                                                temptill -= GlobalMethods.layerthickness_m[row, col, layero];
                                                // GlobalMethods.layerthickness_m[row, col, layero] = 0;
                                                layero++;
                                            }
                                            // Debug.WriteLine("till8");
                                            // transport eroded fraction
                                            frac_eroded = temptill / GlobalMethods.layerthickness_m[row, col, layero];
                                            // mass fraction eroded
                                            for (int tex = 0; tex < 5; tex++)
                                            {
                                                GlobalMethods.texture_kg[row + i, col + j, 0, tex] += GlobalMethods.texture_kg[row, col, layero, tex] * frac_eroded;
                                                GlobalMethods.texture_kg[row, col, layero, tex] -= GlobalMethods.texture_kg[row, col, layero, tex] * frac_eroded;
                                            }
                                            GlobalMethods.young_SOM_kg[row + i, col + j, 0] += GlobalMethods.young_SOM_kg[row, col, layero] * frac_eroded;
                                            GlobalMethods.young_SOM_kg[row, col, layero] -= GlobalMethods.young_SOM_kg[row, col, layero] * frac_eroded;
                                            GlobalMethods.old_SOM_kg[row + i, col + j, 0] += GlobalMethods.old_SOM_kg[row, col, layero] * frac_eroded;
                                            GlobalMethods.old_SOM_kg[row, col, layero] -= GlobalMethods.old_SOM_kg[row, col, layero] * frac_eroded;


                                            GlobalMethods.layerthickness_m[row, col, layero] = thickness_calc(row, col, layero);
                                            GlobalMethods.layerthickness_m[row + i, col + j, 0] = thickness_calc(row, col, layero);

                                            //if necessary, i.e. an entire layer removed, shift cells up
                                            if (layero > 0)
                                            {
                                                try { remove_empty_layers(row, col); update_all_soil_thicknesses(row, col); }
                                                catch { Debug.WriteLine("Error in removing empty layers after tillage"); }
                                            }
                                        }//end if
                                    }//end if novalues
                                }//end if borders
                            }//end for j
                        }//end for i
                    } //end if GlobalMethods.tillfields
                }   // end  for 
                    // Debug.WriteLine("till9");
                    // 3. Update elevation changes
                for (row = 0; row < GlobalMethods.nr; row++)
                {
                    for (col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (GlobalMethods.dtm[row, col] != -9999)
                        {
                            double old_soil_thickness = GlobalMethods.soildepth_m[row, col];
                            update_all_soil_thicknesses(row, col);
                            double new_soil_thickness = total_soil_thickness(row, col);

                            GlobalMethods.dtm[row, col] += new_soil_thickness - old_soil_thickness;
                            GlobalMethods.dtmchange[row, col] += new_soil_thickness - old_soil_thickness;
                            GlobalMethods.sum_tillage[row, col] += new_soil_thickness - old_soil_thickness;
                            GlobalMethods.soildepth_m[row, col] = new_soil_thickness;
                            if (GlobalMethods.till_result[row, col] > 0) { total_sum_tillage += new_soil_thickness - old_soil_thickness; }
                        }
                    }
                }

                total_tillage_statuspanel.Text = string.Format("till {0:F0} * 1000 m3", total_sum_tillage * GlobalMethods.dx * GlobalMethods.dx / 1000);
                // Debug.WriteLine("\n--tillage overview--");
                // Debug.WriteLine(" tilled a total of " + total_sum_tillage * GlobalMethods.dx * GlobalMethods.dx / 1000 + " * 1000 m3");
                double mass_after = total_catchment_mass();
                if (Math.Abs(mass_before - mass_after) > 0.0001)
                {
                    Debug.WriteLine("err_ti3");
                }
            }
            catch
            {
                Debug.WriteLine("err_ti4");

            }
        }

        /*
        private void calculate_creep()
        {
            this.InfoStatusPanel.Text = "GlobalMethods.creep calculation";
            int row, col,
                        i, j,
                        nb_ok;
            double
                        dhmin, dhe_tol, dhs_tol,
                        slope_sum, dhmax, dz_min, GlobalMethods.d_x, dz_max, dh1, dh,
                        fraction,
                        temp, tempcreep, tempdep,
                        slope;


            nb_ok = 0; nb_check = 0; all_grids = 0;
            dhmin = -9999; dhe_tol = 0.000001; dhs_tol = 0.000001;
            for (row = 0; row < GlobalMethods.nr; row++)
            {
                for (col = 0; col < GlobalMethods.nc; col++)
                {
                    GlobalMethods.creep[row, col] = 0;	// neighbour check is 0 is false
                }
            }

            int runner = 0;
            for (runner = GlobalMethods.number_of_data_cells - 1; runner >= 0; runner--)
            {           // the GlobalMethods.index is sorted from low to high values, but flow goes from high to low
                row = row_index[runner]; col = col_index[runner];
                // into loop for surrounding grids of certain grid
                // Start first the slope_sum loop for all lower neighbour grids
                slope_sum = 0; dhmax = 0; dz_min = -9999.99; GlobalMethods.d_x = GlobalMethods.dx;
                dz_max = -1; dhtemp = -99999.99; dhmin = (-9999.99);

                for (i = (-1); i <= 1; i++)
                {
                    for (j = (-1); j <= 1; j++)
                    {
                        dh = 000000; dh1 = 000; dhtemp = -99999.99; GlobalMethods.d_x = GlobalMethods.dx;
                        if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                        {    // boundaries
                            if ((GlobalMethods.dtm)[row + i, col + j] != (-9999))
                            {
                                dh = ((GlobalMethods.dtm)[row, col] - (GlobalMethods.dtm)[row + i, col + j]);
                                if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                                if (dh < 000000)
                                {           // i j is a higher neighbour
                                    if (dh > dz_min) { dz_min = dh; }
                                    if (dh1 > dhmin + dhs_tol) { dhmin = (dh1 + dhs_tol); }
                                }
                                if (dh > 000000)
                                {           // i j is a lower neighbour
                                    if (dh1 > dhmax - dhe_tol) { dhmax = (dh1 - dhe_tol); }
                                    dh = dh / GlobalMethods.d_x;
                                    if (dh > dz_max) { dz_max = dh; }
                                    dh = Math.Pow(dh, conv_fac);
                                    slope_sum = slope_sum + dh;
                                }//end if
                            }//end if novalues
                        }// end if boundaries
                    }//end for
                }//end for
                if (dhmin == -9999.99) { dhmin = 0; } else { dhmin = -dhmin; }
                if (dhmax <= 0.0) { dhmax = 0.0; } else { dhmax = -dhmax; }
                for (i = (-1); i <= 1; i++)
                {
                    for (j = (-1); j <= 1; j++)
                    {
                        dh = 0.000000; fraction = 0.0;
                        frac_dis = 0.0;
                        GlobalMethods.d_x = GlobalMethods.dx;
                        if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0))) // boundaries
                        {
                            if ((GlobalMethods.dtm)[row + i, col + j] != (-9999.0))
                            {
                                dh = ((GlobalMethods.dtm)[row, col] - (GlobalMethods.dtm)[row + i, col + j]);
                                temp = (GlobalMethods.dtm)[row + i, col + j];
                                // Multiple Flow: If there are lower neighbours start evaluating
                                if (dh > 0.000000)
                                {
                                    // fraction of discharge into a neighbour grid
                                    if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                                    slope = dh / GlobalMethods.d_x;
                                    dh = dh / GlobalMethods.d_x;
                                    dh = Math.Pow(dh, conv_fac);
                                    fraction = (dh / slope_sum);
                                    // Tillage erosion calculation
                                    tempcreep = fraction * slope * diffusivity_creep / GlobalMethods.dx;
                                    tempdep = GlobalMethods.soildepth_m[row, col] * -1.0;
                                    if (tempdep > dhmax) tempdep = dhmax;
                                    if (tempcreep < (tempdep)) tempcreep = tempdep;
                                    GlobalMethods.creep[row, col] -= tempcreep;
                                    GlobalMethods.creep[row + i, col + j] += tempcreep;
                                    GlobalMethods.soildepth_m[row, col] -= tempcreep;
                                    GlobalMethods.soildepth_m[row + i, col + j] += tempcreep;
                                    if (GlobalMethods.soildepth_m[row, col] < 0) { GlobalMethods.soildepth_m[row, col] = 0; }
                                    if (GlobalMethods.soildepth_m[row + i, col + j] < 0) { GlobalMethods.soildepth_m[row + i, col + j] = 0; }
                                    GlobalMethods.dtm[row, col] -= tempcreep;
                                    GlobalMethods.dtm[row + i, col + j] += tempcreep;
                                }//end if
                            }//end if novalues
                        }//end if borders
                    }//end for j
                }//end for i
            }		// end for sorted 
            for (row = 0; row < GlobalMethods.nr; row++)
            {
                for (col = 0; col < GlobalMethods.nc; col++)
                {
                    GlobalMethods.dtmchange[row, col] += (GlobalMethods.creep[row, col]);
                    GlobalMethods.sum_creep_grid[row, col] += (GlobalMethods.creep[row, col]);
                    if (GlobalMethods.creep[row, col] > 0) sum_creep += (GlobalMethods.creep[row, col]);
                }
            }
        }
        */ // oude GlobalMethods.creep berekeningen


        /*private void calculate_creep()
        {
            try
            {
                this.InfoStatusPanel.Text = "GlobalMethods.creep calculation";
                int row, col,
                            i, j,
                            nb_ok,
                            NA_dem;
                double
                            dhmin, dhe_tol, dhs_tol,
                            slope_sum, dhmax, dz_min, GlobalMethods.d_x, dz_max, dh1, dh,
                            fraction,
                            temp, tempcreep, tempdep,
                            slope;

                nb_ok = 0; nb_check = 0; all_grids = 0;
                dhmin = -9999; dhe_tol = 0.000001; dhs_tol = 0.000001;
                for (row = 0; row < GlobalMethods.nr; row++)
                {
                    for (col = 0; col < GlobalMethods.nc; col++)
                    {
                        GlobalMethods.creep[row, col] = 0;    // neighbour check is 0 is false
                    }
                }
                NA_dem = NA_in_DEM();
                if (NA_dem != NA_in_DEM()) { Debugger.Break(); }

                int runner = 0;
                for (runner = GlobalMethods.number_of_data_cells - 1; runner >= 0; runner--)
                {           // the GlobalMethods.index is sorted from low to high values, but flow goes from high to low
                    row = row_index[runner]; col = col_index[runner];
                    // into loop for surrounding grids of certain grid
                    // Start first the slope_sum loop for all lower neighbour grids
                    slope_sum = 0; dhmax = 0; dz_min = -9999.99; GlobalMethods.d_x = GlobalMethods.dx;
                    dz_max = -1; dhtemp = -99999.99; dhmin = (-9999.99);
                    // if(row == 31 && col == 12) { Debug.WriteLine("creep1"); displaysoil(row, col); }
                    if (thickness_calc(row, col, 0) < 0) { displaysoil(row, col); Debugger.Break(); }

                    for (i = (-1); i <= 1; i++)
                    {
                        for (j = (-1); j <= 1; j++)
                        {
                            dh = 000000; dh1 = 000; dhtemp = -99999.99; GlobalMethods.d_x = GlobalMethods.dx;
                            if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                            {    // boundaries
                                if ((GlobalMethods.dtm)[row + i, col + j] != (-9999))
                                {
                                    dh = ((GlobalMethods.dtm)[row, col] - (GlobalMethods.dtm)[row + i, col + j]);
                                    if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                                    if (dh < 000000)
                                    {           // i j is a higher neighbour
                                        if (dh > dz_min) { dz_min = dh; }
                                        if (dh1 > dhmin + dhs_tol) { dhmin = (dh1 + dhs_tol); }
                                    }
                                    if (dh > 000000)
                                    {           // i j is a lower neighbour
                                        if (dh1 > dhmax - dhe_tol) { dhmax = (dh1 - dhe_tol); }
                                        dh = dh / GlobalMethods.d_x;
                                        if (dh > dz_max) { dz_max = dh; }
                                        dh = Math.Pow(dh, conv_fac);
                                        slope_sum = slope_sum + dh;
                                    }//end if
                                }//end if novalues
                            }// end if boundaries
                        }//end for j
                    }//end for i
                    if (NA_dem != NA_in_DEM()) { Debugger.Break(); }
                    if (thickness_calc(row, col, 0) < 0) { displaysoil(row, col); Debugger.Break(); }

                    // if (row == 31 && col == 12) { Debug.WriteLine("creep2"); displaysoil(row, col); }

                    if (dhmin == -9999.99) { dhmin = 0; } else { dhmin = -dhmin; }
                    if (dhmax <= 0.0) { dhmax = 0.0; } else { dhmax = -dhmax; }
                    for (i = (-1); i <= 1; i++)
                    {
                        for (j = (-1); j <= 1; j++)
                        {
                            dh = 0.000000; fraction = 0.0;
                            frac_dis = 0.0;
                            GlobalMethods.d_x = GlobalMethods.dx;
                            if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0))) // boundaries
                            {
                                if (NA_dem != NA_in_DEM()) { Debugger.Break(); }
                                if (thickness_calc(row, col, 0) < 0) { displaysoil(row, col); Debugger.Break(); }

                                if ((GlobalMethods.dtm)[row + i, col + j] != (-9999.0))
                                {
                                    dh = ((GlobalMethods.dtm)[row, col] - (GlobalMethods.dtm)[row + i, col + j]);
                                    temp = (GlobalMethods.dtm)[row + i, col + j];
                                    // Multiple Flow: If there are lower neighbours start evaluating
                                    if (dh > 0.000000)
                                    {
                                        // if (row == 31 && col == 12) { Debug.WriteLine("creep3"); displaysoil(row, col); }

                                        // fraction of discharge into a neighbour grid
                                        if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                                        slope = dh / GlobalMethods.d_x;
                                        dh = dh / GlobalMethods.d_x;
                                        dh = Math.Pow(dh, conv_fac);
                                        fraction = (dh / slope_sum);
                                        tempcreep = fraction * slope * diffusivity_creep / GlobalMethods.dx;
                                        tempdep = GlobalMethods.soildepth_m[row, col] * -1.0;
                                        // if (tempcreep > 1) { Debugger.Break(); }
                                        if (tempdep > dhmax) tempdep = dhmax;
                                        if (tempcreep < (tempdep)) tempcreep = tempdep;
                                        if (NA_dem != NA_in_DEM()) { Debugger.Break(); }
                                        if (thickness_calc(row, col, 0) < 0) { displaysoil(row, col); Debugger.Break(); }


                                        //instead, we now call a separate function that calculates for every cell-cell combo the layer-mass of GlobalMethods.creep, stores effects in a big matrix, and then evaluate results of that matrix after the timestep.

                                        //if (row == 31 && col == 12)
                                        //{
                                        //    Debug.WriteLine("creep4"); displaysoil(row, col);
                                        //    Debug.WriteLine("tempcreep = {0}", tempcreep);

                                        //}
                                        if (NA_dem != NA_in_DEM()) { Debugger.Break(); }
                                        // if (thickness_calc(row, col, 0) < 0) { displaysoil(row, col); Debugger.Break(); }

                                        // oldsoildepths
                                        double dsoil_source = total_soil_thickness(row, col);
                                        double dsoil_sink = total_soil_thickness(row + i, col + j);
                                        calc_creep_layers(row, col, i, j, tempcreep);

                                        // if (row == 31 && col == 12) { Debug.WriteLine("creep5"); displaysoil(row, col); }

                                        if (NA_dem != NA_in_DEM()) { Debugger.Break(); }
                                        // if (thickness_calc(row, col, 0) < 0) { displaysoil(row, col); Debugger.Break(); }

                                        // update soil depths
                                        update_all_soil_thicknesses(row, col);
                                        update_all_soil_thicknesses(row + i, col + j);
                                        // map updates

                                        // if (row == 31 && col == 12) { Debug.WriteLine("creep6"); displaysoil(row, col); }

                                        if (NA_dem != NA_in_DEM()) { Debugger.Break(); }
                                        if (thickness_calc(row, col, 0) < 0) { displaysoil(row, col); Debugger.Break(); }

                                        double dz_source = total_soil_thickness(row, col) - dsoil_source; // change in soil depth
                                        double dz_sink = total_soil_thickness(row + i, col + j) - dsoil_sink; // change in soil depth
                                        GlobalMethods.creep[row, col] += dz_source;
                                        GlobalMethods.creep[row + i, col + j] += dz_sink;
                                        GlobalMethods.soildepth_m[row, col] += dz_source;
                                        GlobalMethods.soildepth_m[row + i, col + j] += dz_sink;
                                        if (GlobalMethods.soildepth_m[row, col] < 0) { GlobalMethods.soildepth_m[row, col] = 0; }
                                        if (GlobalMethods.soildepth_m[row + i, col + j] < 0) { GlobalMethods.soildepth_m[row + i, col + j] = 0; }
                                        GlobalMethods.dtm[row, col] += dz_source;
                                        GlobalMethods.dtm[row + i, col + j] += dz_sink;
                                        GlobalMethods.dtmchange[row, col] += dz_source; //MMS
                                        GlobalMethods.dtmchange[row + i, col + j] += dz_sink; //MMS

                                        if (NA_dem != NA_in_DEM()) { Debugger.Break(); }
                                        if (thickness_calc(row, col, 0) < 0) { displaysoil(row, col); Debugger.Break(); }

                                        // if (row == 31 && col == 12) { Debug.WriteLine("creep7"); displaysoil(row, col); }


                                    }//end if
                                }//end if novalues
                            }//end if borders
                        }//end for j
                    }//end for i
                }		// end for sorted 
            }
            catch
            {
                Debugger.Break();
            }

        }
        */ // GlobalMethods.creep compatible with multiple soil layers, but with diffusivity in meters. The new function below calculated movement in kg/m2/y MM

        private void calculate_creep()
        {
            // Debug.WriteLine("start of GlobalMethods.creep");
            try
            {
                if (NA_in_map(GlobalMethods.dtm) > 0 | NA_in_map(GlobalMethods.soildepth_m) > 0)
                {
                    Debug.WriteLine("err_cr1");
                    Debugger.Break();
                }

                guiVariables.InfoStatusPanel = "GlobalMethods.creep calculation";
                int row, col,
                            i, j,
                            nb_ok,
                            NA_dem;
                double
                            dhmin, dhe_tol, dhs_tol,
                            slope_sum, dhmax, dz_min, dz_max, dh1, dh,
                            fraction,
                            temp, tempcreep, tempdep,
                            slope,
                            potential_creep_kg = 0, local_creep_kg = 0;

                nb_ok = 0; nb_check = 0; all_grids = 0;
                dhmin = -9999; dhe_tol = 0.00000; dhs_tol = 0.00000;

                //NA_dem = NA_in_DEM();
                //if (NA_dem != NA_in_DEM()) { Debugger.Break(); }

                int runner = 0;

                for (runner = GlobalMethods.number_of_data_cells - 1; runner >= 0; runner--)
                {           // the GlobalMethods.index is sorted from low to high values, but flow goes from high to low
                    row = GlobalMethods.row_index[runner]; col = GlobalMethods.col_index[runner];
                    if (GlobalMethods.dtm[row, col] != -9999)
                    {
                        //Debug.WriteLine("cr1");
                        // into loop for surrounding grids of certain grid
                        // Start first the slope_sum loop for all lower neighbour grids
                        slope_sum = 0; dhmax = 0; dz_min = -9999.99; GlobalMethods.d_x = GlobalMethods.dx;
                        dz_max = -1; dhtemp = -99999.99; dhmin = (-9999.99);
                        // if(row == 31 && col == 12) { Debug.WriteLine("creep1"); displaysoil(row, col); }
                        if (thickness_calc(row, col, 0) < 0)
                        {
                            displaysoil(row, col); Debug.WriteLine("err_cr2");
                        }

                        for (i = (-1); i <= 1; i++)
                        {
                            for (j = (-1); j <= 1; j++)
                            {
                                dh = 000000; dh1 = 000; dhtemp = -99999.99; GlobalMethods.d_x = GlobalMethods.dx;
                                if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                                {    // boundaries
                                    if ((GlobalMethods.dtm)[row + i, col + j] != (-9999))
                                    {
                                        dh = ((GlobalMethods.dtm)[row, col] - (GlobalMethods.dtm)[row + i, col + j]);
                                        if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                                        if (dh < 000000)
                                        {           // i j is a higher neighbour
                                            if (dh > dz_min) { dz_min = dh; }
                                            if (dh1 > dhmin + dhs_tol) { dhmin = (dh1 + dhs_tol); }
                                        }
                                        if (dh > 000000)
                                        {           // i j is a lower neighbour
                                            if (dh > dhmax - dhe_tol) { dhmax = (dh - dhe_tol); }
                                            dh = dh / GlobalMethods.d_x;
                                            if (dh > dz_max) { dz_max = dh; }
                                            dh = Math.Pow(dh, conv_fac);
                                            slope_sum = slope_sum + dh;
                                        }//end if
                                    }//end if novalues
                                }// end if boundaries
                            }//end for j
                        }//end for i
                         //if (NA_dem != NA_in_DEM()) { Debugger.Break(); }
                        if (thickness_calc(row, col, 0) < 0)
                        {
                            displaysoil(row, col); Debug.WriteLine("err_cr3");
                        }
                        if (dhmax < 0)
                        {
                            Debug.WriteLine("err_cr4");
                        }
                        //Debug.WriteLine("cr2");

                        // calculate potential GlobalMethods.creep in kg
                        double maxslope = Math.Atan(dz_max); // max slope in radians
                        // potential_creep_kg = 4.5;
                        potential_creep_kg = Convert.ToDouble(potential_bioturbation_textbox.Text);

                        if (guiVariables.Daily_water)
                        {
                            if (aridity_vegetation[row, col] < 1) { potential_creep_kg = 4 + 0.3; } // grassland
                            else { potential_creep_kg = 4 + 1.3; } // forest
                                                                   // standard potential GlobalMethods.creep of 4 kg. 0.3 or 1.3 is added, based on vegetation type. Rates are derived from Wilkinson 2009: breaking ground and Gabet
                        }



                        local_creep_kg = potential_creep_kg * Math.Sin(maxslope) * Math.Cos(maxslope) * GlobalMethods.dx * GlobalMethods.dx * dt; //Equation from gabet et al., 2003 https://doi.org/10.1146/annurev.earth.31.100901.141314 
                                                                                                                      //Debug.WriteLine("cr3");

                        if (local_creep_kg > 0)
                        {

                            if (dhmin == -9999.99) { dhmin = 0; } else { dhmin = -dhmin; }
                            if (dhmax <= 0.0) { dhmax = 0.0; } else { dhmax = -dhmax; }
                            for (i = (-1); i <= 1; i++)
                            {
                                for (j = (-1); j <= 1; j++)
                                {
                                    dh = 0.000000; fraction = 0.0;
                                    frac_dis = 0.0;
                                    GlobalMethods.d_x = GlobalMethods.dx;
                                    //if (col == 1 | col == (GlobalMethods.nc - 1))
                                    //{ Debugger.Break(); }
                                    if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0))) // boundaries
                                    {
                                        //if (NA_dem != NA_in_DEM()) { Debugger.Break(); }
                                        if (thickness_calc(row, col, 0) < 0)
                                        {
                                            displaysoil(row, col); Debug.WriteLine("err_cr5");
                                        }

                                        if ((GlobalMethods.dtm)[row + i, col + j] != (-9999.0))
                                        {
                                            dh = (GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + i, col + j]);
                                            temp = GlobalMethods.dtm[row + i, col + j];
                                            // if (NA_dem != NA_in_DEM()) { Debugger.Break(); }
                                            // Multiple Flow: If there are lower neighbours start evaluating
                                            if (dh > 0.000000)
                                            {
                                                // if (row == 31 && col == 12) { Debug.WriteLine("creep3"); displaysoil(row, col); }
                                                //Debug.WriteLine("Cr1, GlobalMethods.dtm = {0}", GlobalMethods.dtm[row, col]);
                                                // fraction of discharge into a neighbour grid
                                                if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                                                slope = dh / GlobalMethods.d_x;
                                                dh = dh / GlobalMethods.d_x;
                                                dh = Math.Pow(dh, conv_fac);
                                                fraction = (dh / slope_sum);
                                                tempcreep = fraction * local_creep_kg; //MM develop. Original function was fraction*slope*diffusivity. Do I need to add slope in calculations?



                                                //// oldsoildepths
                                                double dsoil_source = total_soil_thickness(row, col);
                                                double dsoil_sink = total_soil_thickness(row + i, col + j);
                                                //displaysoil(row + i, col + j);
                                                double oldmass = total_soil_mass(row, col) + total_soil_mass(row + i, col + j);

                                                double oldmass_source = total_soil_mass(row, col);
                                                double oldmass_sink = total_soil_mass(row + i, col + j);

                                                calc_creep_layers(row, col, i, j, tempcreep);


                                                // update soil depths                                               
                                                update_all_soil_thicknesses(row, col);
                                                update_all_soil_thicknesses(row + i, col + j);

                                                //if (NA_dem != NA_in_DEM()) { Debugger.Break(); }
                                                if (thickness_calc(row, col, 0) < 0)
                                                {
                                                    displaysoil(row, col);
                                                    Debug.WriteLine("err_cr6");

                                                    Debug.WriteLine(thickness_calc(row, col, 0));
                                                    ;
                                                }

                                                double dsoil_source_new = total_soil_thickness(row, col);
                                                double dsoil_sink_new = total_soil_thickness(row + i, col + j);

                                                //displaysoil(row + i, col + j);

                                                double dz_source = total_soil_thickness(row, col) - dsoil_source; // change in soil depth
                                                double dz_sink = total_soil_thickness(row + i, col + j) - dsoil_sink; // change in soil depth
                                                double newmass = total_soil_mass(row, col) + total_soil_mass(row + i, col + j);
                                                GlobalMethods.creep[row, col] += dz_source;
                                                GlobalMethods.creep[row + i, col + j] += dz_sink;

                                                GlobalMethods.soildepth_m[row, col] += dz_source;
                                                GlobalMethods.soildepth_m[row + i, col + j] += dz_sink;
                                                // if (GlobalMethods.soildepth_m[row, col] > 5) { Debugger.Break(); }
                                                if (GlobalMethods.soildepth_m[row, col] < 0) { GlobalMethods.soildepth_m[row, col] = 0; }
                                                if (GlobalMethods.soildepth_m[row + i, col + j] < 0) { GlobalMethods.soildepth_m[row + i, col + j] = 0; }

                                                GlobalMethods.dtm[row, col] += dz_source;
                                                GlobalMethods.dtm[row + i, col + j] += dz_sink;

                                                GlobalMethods.dtmchange[row, col] += dz_source; //MMS
                                                GlobalMethods.dtmchange[row + i, col + j] += dz_sink; //MMS

                                                //Debug.WriteLine("Cr5, GlobalMethods.dtm = {0}", GlobalMethods.dtm[row, col]);
                                                //displaysoil(row, col);

                                                if (Math.Abs(oldmass - newmass) > 0.00001)
                                                {
                                                    Debug.WriteLine("err_cr7");

                                                } //MM Qua hoogte lijkt het hieronder nog mis te gaan. De gewichtscheck hierboven gaat wel goed. 
                                                  // Op DTM zijn de verschillen niet te zien, op GlobalMethods.creep[,] wel
                                                if (Math.Abs(dz_source + dz_sink) > 0.001 & GlobalMethods.t > 1)
                                                {
                                                    // Can occur with a thick lower soil layer, due to small changes in depth->bulk density->thickness. Can be 4 cm for a total soil thickness of 100  m.
                                                    Debug.WriteLine("Creep: Thickness erosion and deposition are not approximately equal");
                                                    //displaysoil(row, col); 
                                                    //displaysoil(row + i, col + j); 
                                                    // Debugger.Break();
                                                }
                                                //if (NA_dem != NA_in_DEM()) { Debugger.Break(); }
                                                if (thickness_calc(row, col, 0) < 0)
                                                {
                                                    displaysoil(row, col); Debug.WriteLine("err_cr8");
                                                }

                                                // if (row == 31 && col == 12) { Debug.WriteLine("creep7"); displaysoil(row, col); }


                                            }//end if
                                        }//end if novalues
                                    }//end if borders
                                }//end for j
                            }//end for i
                        } // end potential_creep_kg>0
                    }       // end for sorted 
                } // end runner
                if (NA_in_map(GlobalMethods.dtm) > 0 | NA_in_map(GlobalMethods.soildepth_m) > 0)
                {
                    Debug.WriteLine("err_cr9");
                }

            }
            catch
            {
                Debug.WriteLine("err_cr10");

            }
            // Debug.WriteLine("end of GlobalMethods.creep");
        }

        private void calc_creep_layers(int row1, int col1, int iiii, int jjjj, double mass_export_soil_kg)
        {
            // tempcreep in kg
            try
            {
                //Debug.WriteLine("Cr3.1, GlobalMethods.dtm = {0}", GlobalMethods.dtm[row1, col1]);
                //displaysoil(row1, col1);
                int layerreceiver = 0;
                double creep_depth_decay_constant = Convert.ToDouble(bioturbation_depth_decay_textbox.Text);

                double frac_dz_lay, frac_overlap_lay, upperdepthdonor = 0, lowerdepthdonor = 0, upperdepthreceiver = 0, lowerdepthreceiver = 0, dsoil = 0, upp_z_lay = 0, int_curve_total, int_curve_lay, mass_export_lay_kg;
                bool C_done = false, lastlayer = false;

                dsoil = total_soil_thickness(row1, col1);

                int_curve_total = 1 / (-creep_depth_decay_constant) * Math.Exp(-creep_depth_decay_constant * 0) - 1 / (-creep_depth_decay_constant) * Math.Exp(-creep_depth_decay_constant * dsoil); // integral over depth decay function, from depth 0 to total soil depth
                upperdepthdonor = 0; //  GlobalMethods.dtm[row1, col1]; using 0 leads to a continuous landscapes, instead of a step-wise pattern
                upperdepthreceiver = 0; // GlobalMethods.dtm[row1 + iiii, col1 + jjjj];
                lowerdepthreceiver = upperdepthreceiver - GlobalMethods.layerthickness_m[row1 + iiii, col1 + jjjj, layerreceiver];

                //if (row1 == 0 & col1 == 0)
                //{
                //    displaysoil(row1, col1);
                //    displaysoil(row1 + iiii, col1 + jjjj);

                //}

                //if (row == 31 && col == 12) { Debug.WriteLine("creep4a. tempcreep = {0}",tempcreep); displaysoil(row, col); }

                for (int lay = 0; lay < GlobalMethods.max_soil_layers; lay++) // test per layer where material moves to
                {

                    if (GlobalMethods.layerthickness_m[row1, col1, lay] > 0)
                    {
                        int_curve_lay = 1 / (-creep_depth_decay_constant) * Math.Exp(-creep_depth_decay_constant * upp_z_lay) - 1 / (-creep_depth_decay_constant) * Math.Exp(-creep_depth_decay_constant * (upp_z_lay + GlobalMethods.layerthickness_m[row1, col1, lay]));//integral over depth decay function, from top of layer to bottom of layer
                        upp_z_lay += GlobalMethods.layerthickness_m[row1, col1, lay];
                        lowerdepthdonor = upperdepthdonor - GlobalMethods.layerthickness_m[row1, col1, lay]; // elevation range donor layer  
                        mass_export_lay_kg = mass_export_soil_kg * (int_curve_lay / int_curve_total); // mass to be removed from layer in kg 

                        //frac_dz_lay = (tempcreep * int_curve_lay / int_curve_total) / GlobalMethods.layerthickness_m[row1, col1, lay]; // fraction that has to be removed
                        frac_overlap_lay = 0; // this fraction will be used to correct for partally overlapping layers 

                        //five options: 
                        //              A donor layer is located completely above receiving layer, exchange with air above receiving layer 0,
                        //              B donor layer partially sticks above upper receiving layer, exchange with air above receiving layer 0,
                        // option A and B will not be used, as we consider the transitio between cells as a continuous curve, i.e. not step-wise pattern
                        // -----------------------------------------------------------
                        //              C (partial) overlap with receiving layer higher than donor layer, 
                        //              D receiving layer fully overlapped by (thicker) donor layer,
                        //              E (partial) overlap with receiving layer lower than donor layer,
                        //              F donor layer fully overlapped by (thicker) receiving layer,
                        //              (G exact overlap (which is like B or E, therefore not explicitly modeled)

                        // Options A and B are about donor layers above the surface of the receiving cell.
                        // Options C-F are about subsurface overlaps of layers, working from higher to lower receiving layers.
                        // This enables update of the receiving layer, when the overlap no longer exists
                        // Exchange of mass is done immediately. After this loop layer thicknesses and DTM are updated for the next calculation

                        // OPTION A: donor layer lies completely above receiving layer
                        // Not possible with the curret setup, where upper depth of donor and receiver both are zero, as is the case in continuous landscapes
                        if (lowerdepthdonor >= upperdepthreceiver && layerreceiver == 0)
                        {
                            frac_overlap_lay = 1;
                            creep_transport(row1, col1, lay, row1 + iiii, col1 + jjjj, layerreceiver, mass_export_lay_kg, frac_overlap_lay);
                            // if(row1==0&&col1==0){ Debug.WriteLine("A, layer " +lay+": " + frac_dz_lay * frac_overlap_lay); }
                            // no need to update rieceiving layer number
                        }

                        // OPTION B. donor layer partly rises above surface source layer. exchange with air above receiving layer 0
                        // Not possible with the curret setup, where upper depth of donor and receiver both are zero, as is the case in continuous landscapes

                        if (upperdepthdonor > upperdepthreceiver && lowerdepthdonor < upperdepthreceiver && layerreceiver == 0)
                        {
                            frac_overlap_lay = (upperdepthdonor - upperdepthreceiver) / GlobalMethods.layerthickness_m[row1, col1, lay];
                            creep_transport(row1, col1, lay, row1 + iiii, col1 + jjjj, layerreceiver, mass_export_lay_kg, frac_overlap_lay);
                            // no need to update receiving layer number, because we only look at air exchange. subsurface exchange will be treated later
                            // if (row1 == 0 && col1 == 0) { Debug.WriteLine("B, layer " + lay + ": " + frac_dz_lay * frac_overlap_lay); }

                        }

                        // OPTION C: (partial) overlap with receiving layer located higher than donor layer
                        if (upperdepthdonor <= upperdepthreceiver && lowerdepthdonor <= lowerdepthreceiver && upperdepthdonor > lowerdepthreceiver)
                        {
                            frac_overlap_lay = (upperdepthdonor - lowerdepthreceiver) / GlobalMethods.layerthickness_m[row1, col1, lay];
                            creep_transport(row1, col1, lay, row1 + iiii, col1 + jjjj, layerreceiver, mass_export_lay_kg, frac_overlap_lay);

                            C_done = true;

                            if (lowerdepthdonor <= lowerdepthreceiver && layerreceiver < (GlobalMethods.max_soil_layers - 1)) // update receiving layer to a lower one. only occurs when lowerdepthdonor == lowerdepthreceiver
                            {
                                layerreceiver += 1;
                                upperdepthreceiver = lowerdepthreceiver;
                                lowerdepthreceiver -= GlobalMethods.layerthickness_m[row1 + iiii, col1 + jjjj, layerreceiver];
                            }
                            // if (row1 == 0 && col1 == 0) { Debug.WriteLine("C, layer " + lay + ": " + frac_dz_lay * frac_overlap_lay); }

                        }

                        // OPTION D: receiving layer completely overlapped by (thicker) donor layer
                        while (upperdepthdonor > upperdepthreceiver && lowerdepthdonor < lowerdepthreceiver && lastlayer == false) // while loop, this can occur several times, when the donor layer completely overlaps receiving layers
                        {
                            frac_overlap_lay = (upperdepthreceiver - lowerdepthreceiver) / GlobalMethods.layerthickness_m[row1, col1, lay];
                            creep_transport(row1, col1, lay, row1 + iiii, col1 + jjjj, layerreceiver, mass_export_lay_kg, frac_overlap_lay);
                            // update receiving layer. the next layer can also be overlapped completely by donor layer
                            if (layerreceiver < (GlobalMethods.max_soil_layers - 1))
                            {
                                layerreceiver += 1;
                                upperdepthreceiver = lowerdepthreceiver;
                                lowerdepthreceiver -= GlobalMethods.layerthickness_m[row1 + iiii, col1 + jjjj, layerreceiver];
                            }
                            else
                            {
                                lastlayer = true;
                            }
                            // if (row1 == 0 && col1 == 0) { Debug.WriteLine("D, layer " + lay + ": " + frac_dz_lay * frac_overlap_lay); }

                        }
                        //OPTION E: overlap with receiving layer lower than donor layer  (take care that this does not evaluate to TRUE when C is also TRUE)
                        if (upperdepthdonor >= upperdepthreceiver && lowerdepthdonor >= lowerdepthreceiver && lowerdepthdonor < upperdepthreceiver && C_done == false)
                        {
                            frac_overlap_lay = (upperdepthreceiver - lowerdepthdonor) / GlobalMethods.layerthickness_m[row1, col1, lay];
                            creep_transport(row1, col1, lay, row1 + iiii, col1 + jjjj, layerreceiver, mass_export_lay_kg, frac_overlap_lay);

                            if (lowerdepthdonor <= lowerdepthreceiver && layerreceiver < (GlobalMethods.max_soil_layers - 1)) // update receiving layer to a lower one
                            {
                                layerreceiver += 1;
                                upperdepthreceiver = lowerdepthreceiver;
                                lowerdepthreceiver -= GlobalMethods.layerthickness_m[row1 + iiii, col1 + jjjj, layerreceiver];
                            }
                            //if (row1 == 0 && col1 == 0)
                            //{
                            //    Debug.WriteLine("E, layer " + lay + ": " + frac_dz_lay * frac_overlap_lay);
                            //    Debug.WriteLine("Dupper = {0}, Dlower = {1}, Rupper = {2}, Rlower = {3}", upperdepthdonor, lowerdepthdonor, upperdepthreceiver, lowerdepthreceiver);
                            //}
                        }

                        //OPTION F, donor layer completely overlapped by (thicker) receiver layer
                        if (upperdepthdonor < upperdepthreceiver && lowerdepthdonor > lowerdepthreceiver)
                        {
                            frac_overlap_lay = 1;
                            creep_transport(row1, col1, lay, row1 + iiii, col1 + jjjj, layerreceiver, mass_export_lay_kg, frac_overlap_lay);
                            // no update of receiving layer required
                            // if (row1 == 0 && col1 == 0) { Debug.WriteLine("F, layer " + lay + ": " + frac_dz_lay * frac_overlap_lay); }
                        }

                        //OPTION H, donor soil might be absent. Material moves to upper layer of receiving layer, if elevation allows
                        if (total_soil_thickness(row1 + iiii, col1 + jjjj) == 0) // if receiving cell is bare rock 
                        {
                            bool partial_overlap = true;
                            if ((GlobalMethods.dtm[row1, col1] + upperdepthdonor) < GlobalMethods.dtm[row1 + iiii, col1 + jjjj]) { frac_overlap_lay = 0; partial_overlap = false; } // donor layer lies completely below surface of receiving cell
                            if ((GlobalMethods.dtm[row1, col1] + lowerdepthdonor) > GlobalMethods.dtm[row1 + iiii, col1 + jjjj]) { frac_overlap_lay = 1; partial_overlap = false; } // donor layer lies completely above surface of receiving cell

                            if (partial_overlap == true)
                            {
                                frac_overlap_lay = ((GlobalMethods.dtm[row1, col1] + upperdepthdonor) - GlobalMethods.dtm[row1 + iiii, col1 + jjjj]) / (upperdepthdonor - lowerdepthdonor);
                            } // donor layer lies partially above surface of receiving cell
                            if (frac_overlap_lay > 1) { frac_overlap_lay = 1; } // fraction can be a bit higher due to rounding errors

                            if (Double.IsInfinity(frac_overlap_lay) | Double.IsNaN(frac_overlap_lay) | frac_overlap_lay > 1)
                            {
                                Debug.WriteLine("err_cr11");
                            } // something went wrong in calculating overlapping fraction. Either divided by zero, a non-real answer or a fraction larger than 1


                            if (frac_overlap_lay > 0) // If the donor layer is (partially) above the bare bedrock of the receiving cell, everything can move to next cell:
                            {
                                frac_overlap_lay = 1;
                                creep_transport(row1, col1, lay, row1 + iiii, col1 + jjjj, 0, mass_export_lay_kg, frac_overlap_lay);
                            }
                        }









                        upperdepthdonor = lowerdepthdonor;
                        C_done = false;
                    } // end layerthickness > 0
                } // end layers
                  // if (row == 31 && col == 12) { Debug.WriteLine("creep4b"); displaysoil(row, col); }

            } // end of try
            catch
            {
                Debug.WriteLine("Error in time {0}, row {1}, col{2}, receiving row {4}, col {5}", GlobalMethods.t, row1, col1, row1 + iiii, col1 + jjjj);
                Debug.WriteLine("err_cr12");

            }
        } // end calc_creep_layers

        private void creep_transport(int fromrow, int fromcol, int fromlay, int torow, int tocol, int tolay, double mass_export, double fraction_overlap)
        {
            try
            {
                //Debug.WriteLine("Cr3.1.1, GlobalMethods.dtm = {0}", GlobalMethods.dtm[fromrow, fromcol]);
                //displaysoil(fromrow, fromcol);

                // double fraction_transport = fraction_dz * fraction_overlap;
                double fraction_transport = mass_export / total_layer_mass(fromrow, fromcol, fromlay); // fraction of mass to be exported
                if (fraction_transport > 1) { fraction_transport = 1; }
                if (fraction_transport < 0)
                {
                    fraction_transport = 0; Debug.WriteLine("err_cr13");

                }
                if (fraction_overlap > 1) { fraction_overlap = 1; }
                for (int tex = 0; tex < 5; tex++)
                {
                    GlobalMethods.texture_kg[torow, tocol, tolay, tex] += GlobalMethods.texture_kg[fromrow, fromcol, fromlay, tex] * fraction_transport;
                    GlobalMethods.texture_kg[fromrow, fromcol, fromlay, tex] -= GlobalMethods.texture_kg[fromrow, fromcol, fromlay, tex] * fraction_transport;
                }
                GlobalMethods.young_SOM_kg[torow, tocol, tolay] += GlobalMethods.young_SOM_kg[fromrow, fromcol, fromlay] * fraction_transport;
                GlobalMethods.young_SOM_kg[fromrow, fromcol, fromlay] -= GlobalMethods.young_SOM_kg[fromrow, fromcol, fromlay] * fraction_transport;
                GlobalMethods.old_SOM_kg[torow, tocol, tolay] += GlobalMethods.old_SOM_kg[fromrow, fromcol, fromlay] * fraction_transport;
                GlobalMethods.old_SOM_kg[fromrow, fromcol, fromlay] -= GlobalMethods.old_SOM_kg[fromrow, fromcol, fromlay] * fraction_transport;

                //Debug.WriteLine("Cr3.1.2, GlobalMethods.dtm = {0}", GlobalMethods.dtm[fromrow, fromcol]);
                //displaysoil(fromrow, fromcol);

            }
            catch
            {
                Debug.WriteLine("crashed during GlobalMethods.creep transport calculations");
                Debug.WriteLine("err_cr14");

            }


        }

        private void calculate_tree_fall()
        {
            double tf_mass_before = total_catchment_mass();

            try
            {
                this.InfoStatusPanel.Text = "tree fall calculation";

                bool fallen = false;
                int i_tf = 0, j_tf = 0;
                // double mass_before_tf = total_catchment_mass();
                double exported_mass_tf = 0, old_soil_depth_m = 0, tree_fall_frac_sum, tf_frac_dx;

                double[] tree_fall_mass, tree_fall_om;
                double[,] tree_fall_frac;
                Random rand = new Random(GlobalMethods.t); // GlobalMethods.t as a random seed
                Random falldirection = new Random(GlobalMethods.t);
                Random age_of_trees = new Random(GlobalMethods.t);
                int P_fall = Convert.ToInt32(Math.Round((1 / tf_frequency) / (GlobalMethods.dx * GlobalMethods.dx)));
                // int P_fall = Convert.ToInt32(Math.Round(1730 / GlobalMethods.dx / GlobalMethods.dx)); // 1/P_fall is the chance of tree fall, per m2, that's why we correct for cell size 
                // Debug.WriteLine("elevation of row 57 and col 40 at GlobalMethods.t {0} is {1}", GlobalMethods.t, GlobalMethods.dtm[57, 40]);
                int rowsource = 0, colsource = 0, rowsink = 0, colsink = 0;
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (GlobalMethods.dtm[row, col] != -9999 & aridity_vegetation[row, col] > 1) // if cell exists, and if there is no grass growing
                        {
                            int chance = rand.Next(0, P_fall);
                            // Debug.WriteLine("tf2");
                            //if (row == 31 & col == 31 & GlobalMethods.t == 226) { displaysoil(31, 31); Debugger.Break(); }
                            if (chance == 1) // if a tree falls
                            {
                                rowsource = row;
                                colsource = col;
                                fallen = true;
                                GlobalMethods.treefall_count[row, col] += 1;
                                // if (row == 0 & col == 5) { Debug.WriteLine("tf on GlobalMethods.t {0}", GlobalMethods.t); }
                                int falldir = falldirection.Next(1, 9); // for now a random fall direction. This can be changed as a function of e.g. slope, GlobalMethods.aspect and dominant wind direction
                                                                        // It appears that these factors don'GlobalMethods.t have a dominant effect:https://doi.org/10.3159/10-RA-011.1 
                                                                        // trees can now fall in 8 directions, to all neighbouring cells. Depending on the distance to these cells, sediments will be redistributed.
                                                                        // neighbours:
                                                                        // 1  2  3
                                                                        // 4  X  5
                                                                        // 6  7  8
                                                                        // DEVELOP change surface roughness as function of tree fall, to promote more infiltration?
                                                                        //determine row direction i_tf
                                if (falldir < 4) { i_tf = -1; }
                                else if (falldir < 6) { i_tf = 0; }
                                else if (falldir < 9) { i_tf = 1; }
                                else { MessageBox.Show("error in tree fall. Fall direction Y not known"); }
                                // Debug.WriteLine("tf3");

                                // determine col direction j_tf
                                if (falldir == 1 | falldir == 4 | falldir == 6) { j_tf = -1; }
                                else if (falldir == 2 | falldir == 7) { j_tf = 0; }
                                else if (falldir == 3 | falldir == 5 | falldir == 8) { j_tf = 1; }
                                else { MessageBox.Show("error in tree fall. Fall direction X not known"); }

                                GlobalMethods.d_x = GlobalMethods.dx; if (falldir == 1 | falldir == 3 | falldir == 6 | falldir == 8) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } // determine lateral fall distance

                                if ((row + i_tf) >= 0 & (row + i_tf) < GlobalMethods.nr & (col + j_tf) >= 0 & (col + j_tf) < GlobalMethods.nc)
                                {
                                    if (GlobalMethods.dtm[row + i_tf, col + j_tf] != -9999)
                                    {
                                        dh = (GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + i_tf, col + j_tf]) / GlobalMethods.d_x; // if receiving cell is in catchment, calculate slope
                                    }
                                    else
                                    {
                                        dh = 0;
                                    }
                                }
                                else
                                {
                                    dh = 0; // receiving cell outside catchment, so we assume a slope of 0 percent in order to do the calculations below
                                }

                                double dh_deg = Math.Atan(dh);

                                // growth of tree root system. A spherical growth model is assumed. 
                                // Maximum root system width W_m = 4 m in a circle. This will later be converted to a square surface with the same area
                                // Maximum root system depth D_m = 0.7 m. // From paper Finke on tree fall
                                // These depths are reached after 150 years (typical life span)
                                if (thickness_calc(row, col, 0) < 0)
                                {
                                    Debug.WriteLine("err_tf1");
                                }
                                double W_m, D_m;
                                int tree_age = age_of_trees.Next(0, age_a_max); // selects age between 0 and maximum age of tree, for growth model
                                if (tree_age <= growth_a_max)
                                {
                                    W_m = W_m_max * (3 / 2 * tree_age / growth_a_max - 0.5 * Math.Pow((tree_age / growth_a_max), 3));
                                    D_m = D_m_max * (3 / 2 * tree_age / growth_a_max - 0.5 * Math.Pow((tree_age / growth_a_max), 3));
                                }
                                else
                                {
                                    W_m = W_m_max;
                                    D_m = D_m_max;
                                } // growth formula of trees, giving variable root sizes

                                // Convert spherical surface area to square surface area, to facilitate calculations
                                W_m = Math.Sqrt(Math.PI * Math.Pow(W_m / 2, 2));
                                // Debug.WriteLine("tf4");

                                // Calculation of transported mass over distance
                                int n_affected_cells = 1; // number of affected cells (in a square, so n_cells * n_cells
                                while (W_m > n_affected_cells * GlobalMethods.dx)
                                {
                                    n_affected_cells += 2;
                                }
                                tree_fall_mass = new double[5];
                                tree_fall_om = new double[2];

                                tree_fall_frac = new double[n_affected_cells, n_affected_cells]; // keep track of fractions that are removed from all cells, so that those fractions can be redistributed in teh same pattern one or a few cells outward. 
                                tree_fall_frac_sum = 0;

                                // Debug.WriteLine("N affected cells = {0}", n_affected_cells);

                                // minimaps(row, col); // Base situation

                                for (int ii = (n_affected_cells - 1) / -2; ii <= (n_affected_cells - 1) / 2; ii++)
                                {
                                    for (int jj = (n_affected_cells - 1) / -2; jj <= (n_affected_cells - 1) / 2; jj++)
                                    {
                                        // Debug.WriteLine("tf10");

                                        if (Math.Abs(ii) == (n_affected_cells - 1) / 2 | Math.Abs(jj) == (n_affected_cells - 1) / 2) // cell on the side
                                        {
                                            if (Math.Abs(ii) == (n_affected_cells - 1) / 2 & Math.Abs(jj) == (n_affected_cells - 1) / 2) // all corner cells. Fraction eroded = overlap^2 / GlobalMethods.dx^2
                                            {
                                                tf_frac_dx = Math.Pow(GlobalMethods.dx - ((n_affected_cells * GlobalMethods.dx - W_m) / 2), 2) / Math.Pow(GlobalMethods.dx, 2);
                                            }
                                            else // all other border cells: fraction = overlap * GlobalMethods.dx / GlobalMethods.dx^2
                                            {
                                                tf_frac_dx = ((GlobalMethods.dx - ((n_affected_cells * GlobalMethods.dx - W_m) / 2)) * GlobalMethods.dx) / Math.Pow(GlobalMethods.dx, 2);
                                            }
                                        }
                                        else
                                        {
                                            tf_frac_dx = 1;
                                        }
                                        if (tf_frac_dx > 1) { MessageBox.Show("df_frac_dx > 1"); }
                                        tree_fall_frac[(n_affected_cells - 1) / 2 + ii, (n_affected_cells - 1) / 2 + jj] = tf_frac_dx;
                                        tree_fall_frac_sum += tf_frac_dx;



                                        if (((row + ii) >= 0) && ((row + ii) < GlobalMethods.nr) && ((col + jj) >= 0) && ((col + jj) < GlobalMethods.nc) && GlobalMethods.dtm[row + ii, col + jj] != -9999)
                                        {
                                            // Debug.WriteLine("tf10b");

                                            old_soil_depth_m = GlobalMethods.soildepth_m[row + ii, col + jj];
                                            // if (GlobalMethods.soildepth_m[row + ii, col + jj] < 50) { Debug.WriteLine("tf1 d = {0}, at GlobalMethods.t {1}, r {2}, c {3}", GlobalMethods.soildepth_m[row + ii, col + jj], GlobalMethods.t, row + ii, col + jj); Debugger.Break(); }

                                            // Debug.WriteLine("tf5");

                                            double depth = 0, tf_frac_dz;
                                            int lay = 0;
                                            while (depth < D_m & lay < GlobalMethods.max_soil_layers)
                                            {
                                                // fraction of lowest layer which is incorporated
                                                if (depth + GlobalMethods.layerthickness_m[row + ii, col + jj, lay] <= D_m) { tf_frac_dz = 1; } // fraction of layer taken up by roots, in z direction
                                                else { tf_frac_dz = (D_m - depth) / GlobalMethods.layerthickness_m[row + ii, col + jj, lay]; }

                                                if (tf_frac_dz > 1) { MessageBox.Show("df_frac_dz > 1"); }

                                                // uptake of sediments
                                                for (int tex = 0; tex < 5; tex++)
                                                {
                                                    tree_fall_mass[tex] += GlobalMethods.texture_kg[row + ii, col + jj, lay, tex] * tf_frac_dz * tf_frac_dx;
                                                    GlobalMethods.texture_kg[row + ii, col + jj, lay, tex] -= GlobalMethods.texture_kg[row + ii, col + jj, lay, tex] * tf_frac_dz * tf_frac_dx;
                                                }
                                                tree_fall_om[0] += GlobalMethods.old_SOM_kg[row + ii, col + jj, lay] * tf_frac_dz * tf_frac_dx;
                                                GlobalMethods.old_SOM_kg[row + ii, col + jj, lay] -= GlobalMethods.old_SOM_kg[row + ii, col + jj, lay] * tf_frac_dz * tf_frac_dx;
                                                tree_fall_om[1] += GlobalMethods.young_SOM_kg[row + ii, col + jj, lay] * tf_frac_dz * tf_frac_dx;
                                                GlobalMethods.young_SOM_kg[row + ii, col + jj, lay] -= GlobalMethods.young_SOM_kg[row + ii, col + jj, lay] * tf_frac_dz * tf_frac_dx;
                                                // Debug.WriteLine("tf10c");

                                                // if (thickness_calc(row, col, 0) < 0) { Debugger.Break(); }

                                                // verder
                                                // check if all fractions are calculated correctly
                                                // what if W_m == GlobalMethods.dx, than the fraction is 0; correct by doing 1-fraction? is the eroded fraction still calculate correctly?
                                                // calculations of fraction have to be corrected
                                                // redistribution to a next cell, with the right distance etc 
                                                // solve mass loss tree fall
                                                // GlobalMethods.dtm output has a lot of -9999 instead of no data, check that in R

                                                // update depth and reference layer
                                                depth += GlobalMethods.layerthickness_m[row + ii, col + jj, lay];
                                                lay += 1;
                                            } // end while depth  < Dm
                                              //Debug.WriteLine("Total soil mass: {0}", total_soil_mass(row + ii, col + jj));
                                              //displaysoil(row + ii, col + jj);
                                            remove_empty_layers(row + ii, col + jj);
                                            //Debug.WriteLine("Total soil mass: {0}", total_soil_mass(row + ii, col + jj));
                                            //displaysoil(row + ii, col + jj);

                                            // Debug.WriteLine("tf11");
                                            update_all_soil_thicknesses(row + ii, col + jj);
                                            update_all_soil_thicknesses(row + ii, col + jj); // meij twice, because bulk density depends on depth. This way, the thickness of the empty layers is set to 0 in the first calculation, and used for bulk density in the second calculation

                                            // Elevation change by erosion (removal of material). 
                                            GlobalMethods.soildepth_m[row + ii, col + jj] = total_soil_thickness(row + ii, col + jj);
                                            //if (GlobalMethods.soildepth_m[row + ii, col + jj] < 50) { Debug.WriteLine("tf2 d = {0}, at GlobalMethods.t {1}, r {2}, c {3}", GlobalMethods.soildepth_m[row + ii, col + jj], GlobalMethods.t, row + ii, col + jj); Debugger.Break(); }
                                            if (thickness_calc(row, col, 0) < 0)
                                            {
                                                Debug.WriteLine("err_tf2");
                                            }
                                            // Debug.WriteLine("tf12");

                                            if (old_soil_depth_m - GlobalMethods.soildepth_m[row + ii, col + jj] > 1)
                                            {
                                                Debug.WriteLine("err_tf3");
                                            }

                                            GlobalMethods.dtm[row + ii, col + jj] += GlobalMethods.soildepth_m[row + ii, col + jj] - old_soil_depth_m;
                                            GlobalMethods.dz_treefall[row + ii, col + jj] += GlobalMethods.soildepth_m[row + ii, col + jj] - old_soil_depth_m;
                                            GlobalMethods.dtmchange[row + ii, col + jj] += GlobalMethods.soildepth_m[row + ii, col + jj] - old_soil_depth_m;
                                            // Debug.WriteLine("erosion by tree fall = {0}", GlobalMethods.soildepth_m[row + ii, col + jj] - old_soil_depth_m);
                                            // Debug.WriteLine("tf13");



                                        } // end ii or jj in the area
                                    } // end jj
                                } // end ii
                                  //if (n_affected_cells > 1) { Debugger.Break(); }
                                  // minimaps(row, col); // After erosion


                                // Redistribution of material, deposition
                                // Debug.WriteLine("tf14");
                                double falldist;
                                double dh_rad = dh_deg * (Math.PI / 180);
                                if (dh < 0) // negative slope, tree falls upslope
                                {
                                    falldist = W_m / 2 * (Math.Cos(dh_rad) - Math.Sin(dh_rad)) - D_m / 2 * (Math.Cos(dh_rad) + Math.Sin(dh_rad));
                                }
                                else // positive or zero slope, tree falls downslope
                                {
                                    falldist = W_m / 2 * (Math.Cos(dh_rad) + Math.Sin(dh_rad)) + D_m / 2 * (Math.Sin(dh_rad) - Math.Cos(dh_rad));
                                }

                                // falldist is the distance where the centerpoint of the soil-root mass ends up. The distribution of the material follows the same pattern as the uptake, only shifted a few cells. 
                                // DEVELOP distribute material over different cell, based on dimensions of the soil-root mass
                                // Debug.WriteLine("tf15");
                                int ndist_cells = 0;
                                while (falldist > (ndist_cells + 0.5) * GlobalMethods.dx)
                                {
                                    ndist_cells += 1;
                                }
                                rowsink = row + ndist_cells * i_tf;
                                colsink = col + ndist_cells * j_tf;
                                for (int ii = (n_affected_cells - 1) / -2; ii <= (n_affected_cells - 1) / 2; ii++)
                                {
                                    for (int jj = (n_affected_cells - 1) / -2; jj <= (n_affected_cells - 1) / 2; jj++)
                                    {
                                        // Debug.WriteLine("tf16");

                                        tf_frac_dx = tree_fall_frac[(n_affected_cells - 1) / 2 + ii, (n_affected_cells - 1) / 2 + jj] / tree_fall_frac_sum;
                                        if (((rowsink + ii) >= 0) && ((rowsink + ii) < GlobalMethods.nr) && ((colsink + jj) >= 0) && ((colsink + jj) < GlobalMethods.nc) && GlobalMethods.dtm[rowsink + ii, colsink + jj] != -9999)
                                        {
                                            // Debug.WriteLine("tf17");

                                            for (int tex = 0; tex < 5; tex++)
                                            {
                                                GlobalMethods.texture_kg[rowsink + ii, colsink + jj, 0, tex] += tree_fall_mass[tex] * tf_frac_dx;
                                            }
                                            GlobalMethods.old_SOM_kg[rowsink + ii, colsink + jj, 0] += tree_fall_om[0] * tf_frac_dx;
                                            GlobalMethods.old_SOM_kg[rowsink + ii, colsink + jj, 0] += tree_fall_om[1] * tf_frac_dx;


                                            // elevation change by deposition
                                            old_soil_depth_m = GlobalMethods.soildepth_m[rowsink + ii, colsink + jj];
                                            double ds_1 = GlobalMethods.soildepth_m[rowsink + ii, colsink + jj];
                                            update_all_soil_thicknesses(rowsink + ii, colsink + jj);
                                            double ds_2 = GlobalMethods.soildepth_m[rowsink + ii, colsink + jj];
                                            update_all_soil_thicknesses(rowsink + ii, colsink + jj); // update twice, to approach real BD value, which is depth dependent
                                            double ds_3 = GlobalMethods.soildepth_m[rowsink + ii, colsink + jj];
                                            //if (GlobalMethods.soildepth_m[rowsink + ii, colsink + jj] < 50)
                                            //{
                                            //    Debug.WriteLine("tf3 d = {0}, at GlobalMethods.t {1}, r {2}, c {3}", GlobalMethods.soildepth_m[rowsink + ii, colsink + jj], GlobalMethods.t, rowsink + ii, colsink + jj);
                                            //    Debugger.Break();
                                            //}

                                            GlobalMethods.soildepth_m[rowsink + ii, colsink + jj] = total_soil_thickness(rowsink + ii, colsink + jj);
                                            GlobalMethods.dtm[rowsink + ii, colsink + jj] += GlobalMethods.soildepth_m[rowsink + ii, colsink + jj] - old_soil_depth_m;
                                            GlobalMethods.dz_treefall[rowsink + ii, colsink + jj] += GlobalMethods.soildepth_m[rowsink + ii, colsink + jj] - old_soil_depth_m;
                                            GlobalMethods.dtmchange[rowsink + ii, colsink + jj] += GlobalMethods.soildepth_m[rowsink + ii, colsink + jj] - old_soil_depth_m;
                                            //Debug.WriteLine("deposition by tree fall = {0}", GlobalMethods.soildepth_m[row + ndist_cells * i, col + ndist_cells * j] - old_soil_depth_m);

                                        } // guiVariables.End_time if GlobalMethods.dtm[,] = -9999
                                        else
                                        {
                                            for (int tex = 0; tex < 5; tex++)
                                            {
                                                exported_mass_tf += tree_fall_mass[tex] * tf_frac_dx;
                                            }
                                            exported_mass_tf += tree_fall_om[0] * tf_frac_dx;
                                            exported_mass_tf += tree_fall_om[1] * tf_frac_dx;
                                        }
                                    } // end jj
                                } // end ii
                                  // minimaps(row, col); // After deposition


                            } // end chance ==1 (tree is falling)
                        } // end GlobalMethods.dtm != -9999
                    } // end col
                } // end row
                  // double mass_after_tf = total_catchment_mass() + exported_mass_tf;
                  //if (mass_before_tf != mass_after_tf) { MessageBox.Show("Tree fall mass not equal. difference = "+ (mass_before_tf - mass_after_tf)); }
                if (fallen == true)
                {
                    // Debug.WriteLine("sink cell 1");
                    //displaysoil(rowsink, colsink);

                    //if (GlobalMethods.t == 4)
                    //{
                    //    Debug.WriteLine("Total soil mass: {0}", total_soil_mass(rowsource,colsource));
                    //    displaysoil(rowsource, colsource);
                    //    Debug.WriteLine("Total soil mass: {0}", total_soil_mass(rowsink,colsink));
                    //    displaysoil(rowsink,colsink);
                    //}
                    // Debug.WriteLine("tf1a");
                    for (int r_tf = 0; r_tf < GlobalMethods.nr; r_tf++)
                    {
                        for (int c_tf = 0; c_tf < GlobalMethods.nc; c_tf++)
                        {
                            remove_empty_layers(r_tf, c_tf);
                            remove_empty_layers(r_tf, c_tf);
                            if (total_soil_mass(r_tf, c_tf) <= 0) { Debugger.Break(); }
                            update_all_soil_thicknesses(r_tf, c_tf);
                        }
                    }
                    soil_update_split_and_combine_layers();
                    // Debug.WriteLine("tf1b");

                    //if (GlobalMethods.t == 4)
                    //{
                    //    Debug.WriteLine("Total soil mass: {0}", total_soil_mass(rowsource, colsource));
                    //    displaysoil(rowsource, colsource);
                    //    Debug.WriteLine("Total soil mass: {0}", total_soil_mass(rowsink, colsink));
                    //    displaysoil(rowsink, colsink);
                    //    Debugger.Break();
                    //}
                    if (NA_in_map(GlobalMethods.dtm) > 0)
                    {
                        Debug.WriteLine("err_tf5");
                    }
                    if (NA_in_map(GlobalMethods.soildepth_m) > 0)
                    {
                        Debug.WriteLine("err_tf6");
                    }

                }
                double tf_mass_after = total_catchment_mass() + exported_mass_tf;
                //if (Math.Abs(tf_mass_before - tf_mass_after)>0.001) { Debugger.Break(); }
            }
            catch
            {
                Debug.WriteLine("err_tf7");

            }

        }

        private void calculate_bedrock_weathering()
        {
            // as function of infiltration?
            double Iavg = 0, Imin = 10000000, Imax = 0;
            if (guiVariables.Daily_water)
            {
                int Icount = 0;
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (GlobalMethods.dtm[row, col] != -9999)
                        {
                            if (Imin > Iy[row, col]) { Imin = Iy[row, col]; }
                            if (Imax < Iy[row, col]) { Imax = Iy[row, col]; }
                            Iavg += Iy[row, col];
                            Icount++;
                        }
                    }
                }
                Iavg /= Icount;
            }
            int soil_layer, lowest_soil_layer;
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    if (GlobalMethods.dtm[row, col] != -9999)
                    {

                        // humped
                        if (rockweath_method.SelectedIndex == 0)
                        {
                            GlobalMethods.bedrock_weathering_m[row, col] = P0 * (Math.Exp(-k1 * GlobalMethods.soildepth_m[row, col]) - Math.Exp(-k2 * GlobalMethods.soildepth_m[row, col])) + Pa;

                        }
                        if (rockweath_method.SelectedIndex == 1)
                        {
                            // exponential (Heimsath, Chappell et al., 2000)
                            GlobalMethods.bedrock_weathering_m[row, col] = P0 * (Math.Exp(-k1 * GlobalMethods.soildepth_m[row, col]));
                        }


                        if (rockweath_method.SelectedIndex == 2)
                        {
                            if (guiVariables.Daily_water)
                            {
                                GlobalMethods.bedrock_weathering_m[row, col] = P0 * -k1 * (Iy[row, col] - Imin) / (Imax - Imin);
                            }
                        }


                        GlobalMethods.soildepth_m[row, col] += GlobalMethods.bedrock_weathering_m[row, col]; // this will really be updated at the end of this timestep, but this is a good approximation for the moment
                                                                                 //we also add this amount of coarse material to the lowest layer of our soil
                        soil_layer = 0; lowest_soil_layer = 0;
                        while (GlobalMethods.layerthickness_m[row, col, soil_layer] > 0)
                        {
                            lowest_soil_layer = soil_layer;
                            soil_layer++;
                        }
                        GlobalMethods.texture_kg[row, col, lowest_soil_layer, 0] += GlobalMethods.bedrock_weathering_m[row, col] * 2700 * GlobalMethods.dx * GlobalMethods.dx;   // to go from m (=m3/m2) to kg, we multiply by m2 and by kg/m3
                    }
                }
            }
        }

        private void calculate_tilting()
        {
            guiVariables.InfoStatusPanel = "tilting calculation";
            int row, col;

            for (row = 0; row < GlobalMethods.nr; row++)
            {
                for (col = 0; col < GlobalMethods.nc; col++)
                {
                    if (tilt_location == 0) { GlobalMethods.dtm[row, col] += tilt_intensity * (col / GlobalMethods.nc); }
                    if (tilt_location == 1) { GlobalMethods.dtm[row, col] += tilt_intensity * (row / GlobalMethods.nr); }
                    if (tilt_location == 2) { GlobalMethods.dtm[row, col] += tilt_intensity * ((GlobalMethods.nc - col) / GlobalMethods.nc); }
                    if (tilt_location == 3) { GlobalMethods.dtm[row, col] += tilt_intensity * ((GlobalMethods.nr - row) / GlobalMethods.nr); }
                }
            }
        } //back to the game_clock

        private void calculate_uplift()
        {
            guiVariables.InfoStatusPanel = "uplift calculation";

            int row, col;

            for (row = 0; row < GlobalMethods.nr; row++)
            {
                for (col = 0; col < GlobalMethods.nc; col++)
                {
                    if (lift_location == 0 && row > lift_location) { GlobalMethods.dtm[row, col] += lift_intensity; }
                    if (lift_location == 1 && row > lift_location) { GlobalMethods.dtm[row, col] += lift_intensity; }
                    if (lift_location == 2 && row > lift_location) { GlobalMethods.dtm[row, col] += lift_intensity; }
                    if (lift_location == 3 && row > lift_location) { GlobalMethods.dtm[row, col] += lift_intensity; }
                }
            }
        } //back to the game_clock

        private void calculate_collapse(double max_slope)
        {
            double slope;
            bool last_time_activity = true;
            while (last_time_activity == true)
            //while ( dh > 150;dh / (GlobalMethods.d_x * Math.Sqrt(2)) > 7.5)
            {
                last_time_activity = false;
                for (int row = 0; row < GlobalMethods.nr; row++)
                {        //visit all cells in the DEM and  ...
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        for (i = (-1); i <= 1; i++)
                        {   // maakt een rondje om de cel
                            for (j = (-1); j <= 1; j++)
                            {
                                if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                                {
                                    dh = GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + i, col + j]; // Hoogteverschil					
                                    if ((row != row + i) && (col != col + j)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; } // schuin rakende gridcellen anders Recht rakende gridcellen 				
                                    slope = dh / GlobalMethods.d_x;
                                    dh_tol = GlobalMethods.d_x * max_slope;
                                    if (slope > max_slope)
                                    {
                                        GlobalMethods.dtm[row, col] -= (dh - dh_tol) / 2;
                                        GlobalMethods.dtm[row + i, col + j] += (dh - dh_tol) / 2;
                                        last_time_activity = true;
                                    }
                                } // end height difference                
                            } // end j rondje
                        } // end i rondje
                    } // end col visit
                } // end row visit visit
            } // end while loop
        } // end void collapse ()

        //void update_OSL_age()
        //{
        //    for(row = 0;row<GlobalMethods.nr;row++)
        //    {
        //        for(col = 0;col<GlobalMethods.nc;col++)
        //        {
        //            if(GlobalMethods.dtm[row, col]!=-9999)
        //            {
        //                OSL_age[row,col,0] = 0; // reset surface layer
        //                for(int lay = 1; lay<GlobalMethods.max_soil_layers;lay++)
        //                {
        //                    if(total_layer_mass(row,col,lay)>0)
        //                    {
        //                        OSL_age[row, col, lay] += 1;
        //                    }

        //                }
        //            }
        //        }
        //    }
        //}

        #endregion

        #region Vegetation code

        double[,] aridity_vegetation;

        void determine_vegetation_type()
        {
            aridity_vegetation = new double[GlobalMethods.nr, GlobalMethods.nc];
            double outflow = 0, aridity, averageOF, outflowcells = 0;

            for (int vrow = 0; vrow < GlobalMethods.nr; vrow++)
            {
                for (int vcol = 0; vcol < GlobalMethods.nc; vcol++)
                {

                    if (GlobalMethods.dtm[vrow, vcol] != -9999)
                    {
                        // adjusted Budyko
                        // aridity (water stress) = P/PET. If PET>P, water stress, aridity < 1.
                        // P is replaced by (I+ETa), Incoming water that infiltrates in the cell is captured in I

                        outflow = guiVariables.OFy_m[vrow, vcol, 0] - guiVariables.OFy_m[vrow, vcol, 9];
                        // aridity = (Iy[vrow, vcol] + ETay[vrow, vcol] - outflow) / ET0y[vrow, vcol];
                        aridity = (Iy[vrow, vcol] + ETay[vrow, vcol]) / ET0y[vrow, vcol];

                        // First, we had (I+ET)*(P/(P+O) / PET. But I think the scaling is not necessary. 
                        if (aridity < 0)
                        {
                            Debug.WriteLine("err_vg1");
                        }
                        aridity_vegetation[vrow, vcol] = aridity;


                        if (aridity < 1)
                        {
                            vegetation_type[vrow, vcol] += 1; // arid / grass
                        }
                        else
                        {
                            vegetation_type[vrow, vcol] += 1000; // humid / forest
                        }
                    }
                }
            }
        }


        double[,] veg_correction_factor;
        void change_vegetation_parameters()
        {
            // vegetation coefficients for ET
            for (int vrow = 0; vrow < GlobalMethods.nr; vrow++)
            {

                //open threads 
                for (int vcol = 0; vcol < GlobalMethods.nc; vcol++)
                {
                    if (GlobalMethods.dtm[vrow, vcol] != -9999)
                    {
                        if (aridity_vegetation[vrow, vcol] < 1) { veg_correction_factor[vrow, vcol] = .75; } // all year long, according to FAO report 56
                        else { veg_correction_factor[vrow, vcol] = .85; } // I took the mid-season coefficient (95) of most deciduous crops and decreased it to 85 to account for less vegetation in other times of the year
                        if (GlobalMethods.t >= (guiVariables.End_time - 500)) { veg_correction_factor[vrow, vcol] = .45; }  // if there is agriculture

                    }
                }

                //join threads
            }

            //// bioturbation andd GlobalMethods.creep depth, same parameter
            //if (vegetation_type[vrow, vcol] == 1) { bioturbation_depth_decay_constant = .75; } //
            //if (vegetation_type[vrow, vcol] == 2) { bioturbation_depth_decay_constant = .85; } // I 




        }



        void calculate_TPI(int windowsize)
        {
            try
            {
                //Debug.WriteLine("Started calculating TPI");
                // check if window size is an uneven number, so the window has a center cell
                if (windowsize % 2 == 0) { MessageBox.Show("window size for TPI calculations should be an uneven number"); }

                int windowrange = (windowsize - 1) / 2;

                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (GlobalMethods.dtm[row, col] != -9999)
                        {
                            double tpisum = 0;
                            double tpicount = 0;

                            // calculate moving window average
                            for (int rr = windowrange * -1; rr <= windowrange; rr++)
                            {
                                for (int cc = windowrange * -1; cc <= windowrange; cc++)
                                {
                                    if (row + rr >= 0 & row + rr < GlobalMethods.nr & col + cc >= 0 & col + cc < GlobalMethods.nc) // if cell exists in the DEM, 
                                    {
                                        if (GlobalMethods.dtm[row + rr, col + cc] != -9999 & rr != 0 & cc != 0) // if cell contains a value and cell isn'GlobalMethods.t the target cell, it's considered in the TPI
                                        {
                                            tpisum += GlobalMethods.dtm[row + rr, col + cc];
                                            tpicount += 1;
                                        }
                                    }
                                }
                            }
                            GlobalMethods.tpi[row, col] = GlobalMethods.dtm[row, col] - (tpisum / tpicount);
                        }

                    }
                }
                //Debug.WriteLine("Finished calculating TPI");
            }
            catch
            {
                Debug.WriteLine("Error in calculating TPI");
            }
        }

        #endregion

        #region analysis after simulation code

        void filter_timeseries(double[] series, int window)
        {
            double[] temp_series = new double[series.Length];
            int iterator;
            int one_sided_window = (window - 1) / 2;
            double sum;
            int adder;
            //for the first couple numbers
            for (iterator = 0; iterator < one_sided_window; iterator++)
            {
                sum = 0;
                for (adder = 0; adder <= iterator + one_sided_window; adder++)
                {
                    sum += series[adder];
                }
                temp_series[iterator] = sum / adder;
            }
            //for the (large) middle section of series
            for (iterator = one_sided_window; iterator < (series.Length - one_sided_window); iterator++)
            {
                sum = 0;
                for (adder = iterator - one_sided_window; adder <= iterator + one_sided_window; adder++)
                {
                    sum += series[adder];
                }
                temp_series[iterator] = sum / window;
            }
            //for the last couple datapoints
            for (iterator = (series.Length - one_sided_window); iterator < series.Length; iterator++)
            {
                sum = 0; int numbers = 0;
                for (adder = iterator - one_sided_window; adder <= series.Length; adder++)
                {
                    sum += series[adder]; numbers++;
                }
                temp_series[iterator] = sum / numbers;
            }
            //paste into original series
            for (iterator = 0; iterator < series.Length; iterator++)
            {
                Debug.WriteLine(" replacing " + series[iterator] + " with " + temp_series[iterator]);
                series[iterator] = temp_series[iterator];
            }

        }

        void calc_write_fourier_transform(double[] timeseries_signal)
        {
            filter_timeseries(timeseries_signal, 9);
            Complex[] signal = new Complex[Convert.ToInt32(guiVariables.End_time)];
            // we have made an empty array of complex numbers which we will fill with the information in our input timeseries   
            for (int i = 0; i < guiVariables.End_time; i++) { signal[i] = timeseries_signal[i]; }
            Transform.FourierForward(signal);
            // we have now transformed the raw signal to frequency space, and signal holds that information. Now what?
            for (int i = 0; i < guiVariables.End_time; i++)
            {
                Debug.WriteLine(signal[i]);
            }
        }

        #endregion

        #region experimental and maintenance code

        void minimaps(int row, int col)
        {
            int lowerrow, upperrow, lowercol, uppercol, disrow, discol;
            lowerrow = row - 4; if (lowerrow < 0) { lowerrow = 0; }
            upperrow = row + 4; if (upperrow > GlobalMethods.nr - 1) { upperrow = GlobalMethods.nr - 1; }
            lowercol = col - 4; if (lowercol < 0) { lowercol = 0; }
            uppercol = col + 4; if (uppercol > GlobalMethods.nc - 1) { uppercol = GlobalMethods.nc - 1; }
            string mess;

            // GlobalMethods.dtm
            Debug.Write(" \n"); Debug.Write("      DEM");
            for (discol = lowercol; discol < (uppercol + 1); discol++)
            {
                //Qs = {0:F8}", tomsedi * GlobalMethods.dx * GlobalMethods.dx
                mess = String.Format("  {0:D10}", discol); Debug.Write(mess);
            }
            Debug.Write(" \n"); mess = String.Format(" {0:D8}", lowerrow); Debug.Write(mess);
            for (disrow = lowerrow; disrow < (upperrow + 1); disrow++)
            {
                for (discol = lowercol; discol < (uppercol + 1); discol++)
                {
                    mess = String.Format("  {0:F6}", GlobalMethods.dtm[disrow, discol]); Debug.Write(mess);
                }
                Debug.Write(" \n"); if ((disrow + 1) <= upperrow) { mess = String.Format(" {0:D8}", (disrow + 1)); Debug.Write(mess); }
            }
            /*
            //water flow
            Debug.Write(" \n"); Debug.Write("      Q");
            for (discol = lowercol; discol < (uppercol + 1); discol++)
            {
                //Qs = {0:F8}", tomsedi * GlobalMethods.dx * GlobalMethods.dx
                mess = String.Format("  {0:D10}", discol); Debug.Write(mess);
            }
            Debug.Write(" \n"); mess = String.Format(" {0:D8}", lowerrow); Debug.Write(mess);
            for (disrow = lowerrow; disrow < (upperrow + 1); disrow++)
            {
                for (discol = lowercol; discol < (uppercol + 1); discol++)
                {
                    mess = String.Format("  {0:F6}", guiVariables.OFy_m[disrow, discol,0]); Debug.Write(mess);
                }
                Debug.Write(" \n"); if ((disrow + 1) <= upperrow) { mess = String.Format(" {0:D8}", (disrow + 1)); Debug.Write(mess); }
            }
            */
            Debug.Write(" \n"); Debug.Write("SedI_TRA_kg");
            for (discol = lowercol; discol < (uppercol + 1); discol++)
            {
                mess = String.Format("  {0:D10}", discol); Debug.Write(mess);
            }
            Debug.Write(" \n"); mess = String.Format(" {0:D10}", lowerrow); Debug.Write(mess);
            for (disrow = lowerrow; disrow < (upperrow + 1); disrow++)
            {
                for (discol = lowercol; discol < (uppercol + 1); discol++)
                {
                    mess = String.Format("  {0:F8}", GlobalMethods.sediment_in_transport_kg[disrow, discol, 0]); Debug.Write(mess);
                }
                Debug.Write(" \n"); if ((disrow + 1) <= upperrow) { mess = String.Format(" {0:D10}", (disrow + 1)); Debug.Write(mess); }
            }

            Debug.Write(" \n"); Debug.Write("fillheightA_m");
            for (discol = lowercol; discol < (uppercol + 1); discol++)
            {
                mess = String.Format("  {0:D10}", discol); Debug.Write(mess);
            }
            Debug.Write(" \n"); mess = String.Format(" {0:D10}", lowerrow); Debug.Write(mess);
            for (disrow = lowerrow; disrow < (upperrow + 1); disrow++)
            {
                for (discol = lowercol; discol < (uppercol + 1); discol++)
                {
                    mess = String.Format("  {0:F8}", GlobalMethods.dtmfill_A[disrow, discol]); Debug.Write(mess);
                }
                Debug.Write(" \n"); if ((disrow + 1) <= upperrow) { mess = String.Format(" {0:D10}", (disrow + 1)); Debug.Write(mess); }
            }

            Debug.Write(" \n"); Debug.Write("GlobalMethods.dz_ero_m    ");
            for (discol = lowercol; discol < (uppercol + 1); discol++)
            {
                mess = String.Format("  {0:D10}", discol); Debug.Write(mess);
            }
            Debug.Write(" \n"); mess = String.Format(" {0:D10}", lowerrow); Debug.Write(mess);
            for (disrow = lowerrow; disrow < (upperrow + 1); disrow++)
            {
                for (discol = lowercol; discol < (uppercol + 1); discol++)
                {
                    mess = String.Format("  {0:F8}", GlobalMethods.dz_ero_m[disrow, discol]); Debug.Write(mess);
                }
                Debug.Write(" \n"); if ((disrow + 1) <= upperrow) { mess = String.Format(" {0:D10}", (disrow + 1)); Debug.Write(mess); }
            }

            Debug.Write(" \n"); Debug.Write(" DEPRESSION");
            for (discol = lowercol; discol < (uppercol + 1); discol++)
            {
                mess = String.Format("  {0:D10}", discol); Debug.Write(mess);
            }
            Debug.Write(" \n"); mess = String.Format(" {0:D10}", lowerrow); Debug.Write(mess);
            for (disrow = lowerrow; disrow < (upperrow + 1); disrow++)
            {
                for (discol = lowercol; discol < (uppercol + 1); discol++)
                {
                    mess = String.Format("  {0:D10}", GlobalMethods.depression[disrow, discol]); Debug.Write(mess);
                }
                Debug.Write(" \n"); if ((disrow + 1) <= upperrow) { mess = String.Format(" {0:D10}", (disrow + 1)); Debug.Write(mess); }
            }

            Debug.Write(" \n"); Debug.Write("    status");
            for (discol = lowercol; discol < (uppercol + 1); discol++)
            {
                mess = String.Format("  {0:D10}", discol); Debug.Write(mess);
            }
            Debug.Write(" \n"); mess = String.Format(" {0:D10}", lowerrow); Debug.Write(mess);
            for (disrow = lowerrow; disrow < (upperrow + 1); disrow++)
            {
                for (discol = lowercol; discol < (uppercol + 1); discol++)
                {
                    mess = String.Format("  {0:D10}", GlobalMethods.status_map[disrow, discol]); Debug.Write(mess);
                }
                Debug.Write(" \n"); if ((disrow + 1) <= upperrow) { mess = String.Format(" {0:D10}", (disrow + 1)); Debug.Write(mess); }
            }



        }

        private void display_thick(int row, int col)
        {

        }

        private void displaysoil(int row, int col)
        {

            int layer; double cumthick = 0; double depth = 0, z_layer = GlobalMethods.dtm[row, col];
            //  if (GlobalMethods.t == 0) { Debug.WriteLine("digitally augering and analysing at row " + row + " col " + col); }//header
            Debug.WriteLine("row col GlobalMethods.t nlayer cumth(m)  thick(m)  depth(m) z(m) coarse(kg) sand(kg)   silt(kg)   clay(kg)   fine(kg)   YOM(kg)    OOM(kg)   YOM/OOM   w%coarse   w%sand   w%silt   w%clay   w%fineclay BD"); //header


            for (layer = 0; layer < GlobalMethods.max_soil_layers; layer++) // only the top layer
            {
                //if (GlobalMethods.layerthickness_m[row, col, layer] > 0)
                //{
                cumthick += GlobalMethods.layerthickness_m[row, col, layer];
                depth -= GlobalMethods.layerthickness_m[row, col, layer] / 2;
                double totalweight = GlobalMethods.texture_kg[row, col, layer, 0] + GlobalMethods.texture_kg[row, col, layer, 1] + GlobalMethods.texture_kg[row, col, layer, 2] + GlobalMethods.texture_kg[row, col, layer, 3] + GlobalMethods.texture_kg[row, col, layer, 4] + GlobalMethods.young_SOM_kg[row, col, layer] + GlobalMethods.old_SOM_kg[row, col, layer];
                try { Debug.WriteLine(row + " " + col + " " + GlobalMethods.t + " " + layer + " " + cumthick + " " + GlobalMethods.layerthickness_m[row, col, layer] + " " + depth + " " + z_layer + " " + GlobalMethods.texture_kg[row, col, layer, 0] + " " + GlobalMethods.texture_kg[row, col, layer, 1] + " " + GlobalMethods.texture_kg[row, col, layer, 2] + " " + GlobalMethods.texture_kg[row, col, layer, 3] + " " + GlobalMethods.texture_kg[row, col, layer, 4] + " " + GlobalMethods.young_SOM_kg[row, col, layer] + " " + GlobalMethods.old_SOM_kg[row, col, layer] + " " + GlobalMethods.young_SOM_kg[row, col, layer] / GlobalMethods.old_SOM_kg[row, col, layer] + " " + GlobalMethods.texture_kg[row, col, layer, 0] / totalweight + " " + GlobalMethods.texture_kg[row, col, layer, 1] / totalweight + " " + GlobalMethods.texture_kg[row, col, layer, 2] / totalweight + " " + GlobalMethods.texture_kg[row, col, layer, 3] / totalweight + " " + GlobalMethods.texture_kg[row, col, layer, 4] / totalweight + " " + GlobalMethods.bulkdensity[row, col, layer]); }
                catch { Debug.WriteLine("Cannot write soilprofile"); }
                depth -= GlobalMethods.layerthickness_m[row, col, layer] / 2;
                z_layer -= GlobalMethods.layerthickness_m[row, col, layer];
                //}
            }

            /*if (GlobalMethods.t < guiVariables.End_time )
            {
                
                for (layer = 0; layer < GlobalMethods.max_soil_layers; layer++) // all layers
                {
                    if (GlobalMethods.layerthickness_m[row, col, layer] > 0)
                    {

                        cumthick += GlobalMethods.layerthickness_m[row, col, layer];
                        double totalweight = GlobalMethods.texture_kg[row, col, layer, 0] + GlobalMethods.texture_kg[row, col, layer, 1] + GlobalMethods.texture_kg[row, col, layer, 2] + GlobalMethods.texture_kg[row, col, layer, 3] + GlobalMethods.texture_kg[row, col, layer, 4] + GlobalMethods.young_SOM_kg[row, col, layer] + GlobalMethods.old_SOM_kg[row, col, layer];
                        try { Debug.WriteLine(GlobalMethods.t + " " + cumthick + " " + GlobalMethods.layerthickness_m[row, col, layer] + " " + GlobalMethods.texture_kg[row, col, layer, 0] + " " + GlobalMethods.texture_kg[row, col, layer, 1] + " " + GlobalMethods.texture_kg[row, col, layer, 2] + " " + GlobalMethods.texture_kg[row, col, layer, 3] + " " + GlobalMethods.texture_kg[row, col, layer, 4] + " " + GlobalMethods.young_SOM_kg[row, col, layer] + " " + GlobalMethods.old_SOM_kg[row, col, layer] + " " + GlobalMethods.young_SOM_kg[row, col, layer] / GlobalMethods.old_SOM_kg[row, col, layer] + " " + (GlobalMethods.texture_kg[row, col, layer, 3] + GlobalMethods.texture_kg[row, col, layer, 4]) / totalweight + " " + GlobalMethods.texture_kg[row, col, layer, 2] / totalweight + " " + GlobalMethods.texture_kg[row, col, layer, 1] / totalweight); }
                        catch { Debug.WriteLine("Cannot write soilprofile"); }
                    }
                }
            }*/
            //Debug.WriteLine("");
        }

        bool check_negative_weight(int row, int col)
        {
            bool check = false;
            for (int layer1 = 0; layer1 < GlobalMethods.max_soil_layers; layer1++)
            {
                for (int tex = 0; tex < 5; tex++)
                {
                    if (GlobalMethods.texture_kg[row, col, layer1, tex] < 0)
                    {
                        check = true;
                    }
                }
            }
            return (check);
        }
        double calc_thickness_from_mass(double[] textures_kg, double yom_kg, double oom_kg)
        {
            //pdf goes here
            double thickness_m = 0, soil_mass_kg = 0;
            double sand_fraction, silt_fraction, combined_fraction, bulk_density;
            //first calculate total soil mass to calculate mass percentages for the fine earth fractions (excluding coarse)
            for (int ir = 1; ir < 5; ir++)
            {
                soil_mass_kg += textures_kg[ir];
            }
            soil_mass_kg += oom_kg + yom_kg;
            sand_fraction = textures_kg[1] / soil_mass_kg;
            silt_fraction = textures_kg[2] / soil_mass_kg;

            //calculate bulk density
            bulk_density = bulk_density_calc(textures_kg[0], textures_kg[1], textures_kg[2], textures_kg[3], textures_kg[4], oom_kg, yom_kg, 0.001); // MM depth of 1 micrometer, because a depth of 0 will result in infinite numbers 
            thickness_m = (soil_mass_kg + textures_kg[0]) / (GlobalMethods.dx * GlobalMethods.dx * bulk_density);  // thickness in m per unit area

            return thickness_m;
        }

        #endregion

        #region depression code

        void findsinks()
        {
            this.InfoStatusPanel.Text = "findtrouble has been entered";
            int number, twoequals = 0, threeequals = 0, moreequals = 0;
            int[] intoutlet = new int[9];
            int x;
            numsinks = 0;
            int row, col;

            for (row = 0; row < GlobalMethods.nr; row++)
            {        //visit all cells in the DEM and  ...
                for (col = 0; col < GlobalMethods.nc; col++)
                {
                    if (GlobalMethods.dtm[row, col] != -9999)
                    {
                        dh = 0.0; high = 0; low = 0; equal = 0; GlobalMethods.status_map[row, col] = 0; number = 0;
                        for (x = 0; x < 9; x++) { intoutlet[x] = 99; }
                        for (i = (-1); i <= 1; i++)
                        {       //make a circle around every cell  ...
                            for (j = (-1); j <= 1; j++)
                            {
                                if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)) && GlobalMethods.dtm[row + i, col + j] != -9999)
                                { //boundaries of grid
                                    number++;
                                    dh = GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + i, col + j];
                                    if (dh < 0) { high++; intoutlet[number - 1] = -1; }     // add one to the higher-nbour-counter and add 'h' to the outlet-string
                                    if (dh > 0) { low++; intoutlet[number - 1] = 1; }       // add one to the lower-nbour-counter and add 'l' to the outlet-string
                                    if (dh == 0) { equal++; intoutlet[number - 1] = 0; }      // add one to the equal-nbour-counter and add 'e' to the outlet-string
                                    //Debug.WriteLine("cell %d %d, alt=%.2f nb %d %d, alt=%.2f dh %.3f low %d, equal %d, high %d\n",row,col,GlobalMethods.dtm[row,col],row+i,col+j,GlobalMethods.dtm[row+i,col+j],dh,low,equal,high); 
                                }  //end if within boundaries
                            }  // end for j
                        } // end for i, we have considered the circle around the cell and counted higher, equal and, therefore, lower neighbours

                        if (low == 0 && intoutlet[7] != 99) { GlobalMethods.status_map[row, col] = 1; numsinks++; }       // if 0 lower cells cells are present and 8 cells in total, we have a sink
                        if (low == 8) { GlobalMethods.status_map[row, col] = -1; }       // if 8 lower cells are present, we have a top
                        if (equal == 1) { twoequals++; }        // this case is rare in real DEMs
                        if (equal == 2) { threeequals++; }      // this case is very rare in real DEMs
                        if (equal > 2) { moreequals++; }        // this case is extremely rare in real DEMs
                    } //end for nodata
                    else
                    {
                        GlobalMethods.status_map[row, col] = 3;
                    }
                }   // end for col
            }  // end for row

            /*
            last but not least we give GlobalMethods.status_map a 3 for all cells on the edge of the DEM, so we can end formation of depressions there
            */

            for (col = 0; col < GlobalMethods.nc; col++)
            {
                GlobalMethods.status_map[0, col] = 3; GlobalMethods.status_map[GlobalMethods.nr - 1, col] = 3;
            }
            for (row = 0; row < GlobalMethods.nr; row++)
            {
                GlobalMethods.status_map[row, 0] = 3; GlobalMethods.status_map[row, GlobalMethods.nc - 1] = 3;
            }


            //reports

            this.InfoStatusPanel.Text = "found " + numsinks + " true sinks in " + GlobalMethods.nr * GlobalMethods.nc + "  cells";
            Debug.WriteLine("\n\n--sinks overview at GlobalMethods.t = " + GlobalMethods.t + "--");

            if (numsinks / (GlobalMethods.nr * GlobalMethods.nc) > 0.0075) { Debug.WriteLine("this DEM contains " + numsinks + " true sinks in " + GlobalMethods.nr * GlobalMethods.nc + "  cells\n That's a lot!"); }
            else { Debug.WriteLine("GlobalMethods.t" + GlobalMethods.t + " this DEM contains " + numsinks + " true sinks in " + GlobalMethods.nr * GlobalMethods.nc + "  cells"); }
            //Debug.WriteLine(" equals: " + twoequals / 2 + " double, " + threeequals / 3 + " triple and about " + moreequals + " larger\n");

            // 
            //.WriteLine(" numberofweirdoutlets: " + numberofweirdoutlets + "\n\n"); 
            // 

        }

        void searchdepressions()
        {
            int z;
            this.InfoStatusPanel.Text = "searchdepressions has been entered";
            for (row = 0; row < GlobalMethods.nr; row++)
            {        //visit all cells in the DEM and  ...
                for (col = 0; col < GlobalMethods.nc; col++)
                {
                    depression[row, col] = 0;     // set depression to zero
                }
            }

            for (int z = 0; z < GlobalMethods.numberofsinks; z++)
            {         // the maximum number of depressions is the number of sinks
                GlobalMethods.drainingoutlet_col[z, 0] = -1;
                GlobalMethods.drainingoutlet_col[z, 0] = -1;
                GlobalMethods.drainingoutlet_col[z, 1] = -1;
                GlobalMethods.drainingoutlet_col[z, 1] = -1;
                GlobalMethods.drainingoutlet_col[z, 2] = -1;
                GlobalMethods.drainingoutlet_col[z, 2] = -1;
                GlobalMethods.drainingoutlet_col[z, 3] = -1;
                GlobalMethods.drainingoutlet_col[z, 3] = -1;
                GlobalMethods.drainingoutlet_col[z, 4] = -1;
                GlobalMethods.drainingoutlet_col[z, 4] = -1;
                depressionlevel[z] = 0;
                depressionsize[z] = 0;
                depressionvolume_m[z] = 0;
                iloedge[z] = 0;
                jloedge[z] = 0;
                iupedge[z] = 0;
                jupedge[z] = 0;
            }

            totaldepressions = 0; totaldepressionsize = 0; maxsize = 0; totaldepressionvolume = 0; largestdepression = -1;
            depressionnumber = 0;
            for (row = 0; row < GlobalMethods.nr; row++)
            {        //visit all cells in the DEM and  ...
                for (col = 0; col < GlobalMethods.nc; col++)
                {
                    if (GlobalMethods.status_map[row, col] == 1 && depression[row, col] == 0)
                    {   // sink  -NODATA cells are never sinks, no need to exclude them explicitly here
                        numberoflowestneighbours = 0;
                        for (lowestneighbourcounter = 0; lowestneighbourcounter < maxlowestnbs; lowestneighbourcounter++)
                        {
                            rowlowestnb[lowestneighbourcounter] = -1;
                            collowestnb[lowestneighbourcounter] = -1;
                        }
                        depressionnumber++;                 // so depressionnumber 0 is not used, neither is depression[r,c] = 0
                        if (depressionnumber == 1153000) { diagnostic_mode = 1; } else { diagnostic_mode = 0; }
                        //Debug.WriteLine(" depressionvolume of depression " + depressionnumber + " is initially " + depressionvolume[depressionnumber]); 
                        totaldepressions++;
                        depressionlevel[depressionnumber] = GlobalMethods.dtm[row, col];
                        iloedge[depressionnumber] = row - 1;
                        iupedge[depressionnumber] = row + 1;
                        jloedge[depressionnumber] = col - 1;
                        jupedge[depressionnumber] = col + 1;
                        if (diagnostic_mode == 1)
                        {
                            Debug.WriteLine(" Sink " + depressionnumber + " located: " + row + "," + col + " alt " + GlobalMethods.dtm[row, col]);
                            Debug.WriteLine(" edges : " + iloedge[depressionnumber] + ", " + iupedge[depressionnumber] + ", " + jloedge[depressionnumber] + ", " + jupedge[depressionnumber]);
                        }
                        depressionsize[depressionnumber] = 1;
                        depression[row, col] = depressionnumber;
                        if (depression[row, col] < 0)
                        {
                            MessageBox.Show("Depression error at row " + row + " col " + col + " dep " + depression[row, col]);
                        }
                        iupradius = 1; jupradius = 1; iloradius = 1; jloradius = 1;
                        depressionready = 0; depressiondrainsout = 0;
                        while (depressionready != 1)
                        {
                            if (depressionnumber == 1153000) { diagnostic_mode = 1; }
                            minaltidiff = -99999999; int already_lower_than_lakelevel = 0;
                            for (i = (-1 * iloradius); i <= iupradius; i++)
                            {       //make a circle around the current cell that is so large that it covers all neighbours of all cells currently in the depression
                                for (j = (-1 * jloradius); j <= jupradius; j++)
                                {
                                    if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0))) // in this case, we DO want to find cells with NODATA, because these are valid end-points of lakes
                                    {
                                        //if (diagnostic_mode == 1) { Debug.WriteLine("now at " + (row + i) + ", " + (col + j)); }
                                        nbismemberofdepression = 0;
                                        for (alpha = (-1); alpha <= 1; alpha++)
                                        {       // check in circle around the current cell if one of its neighbours is member of the depression
                                            for (beta = (-1); beta <= 1; beta++)
                                            {
                                                if (((row + i + alpha) >= 0) && ((row + i + alpha) < GlobalMethods.nr) && ((col + j + beta) >= 0) && ((col + j + beta) < GlobalMethods.nc) && !((alpha == 0) && (beta == 0)))
                                                {
                                                    if (depression[row + i + alpha, col + j + beta] == depressionnumber)
                                                    {   //and only if the cell is an actual neighbour of a cell that is member of the depression
                                                        nbismemberofdepression = 1;
                                                    } // end if nb = member of depression
                                                } // end if boundary
                                            }
                                        } // end for first alpha-beta circle to see if any nb = member of depression. It could have been no neighbour of the current depression but still within ilorarius, iupradius etc.
                                        if (nbismemberofdepression == 1 && depression[row + i, col + j] != depressionnumber)
                                        {   // only in case cell row+i, col+j  has a nb that is already member, we are interested in its altitude diff with depressionlevel
                                            altidiff = depressionlevel[depressionnumber] - GlobalMethods.dtm[row + i, col + j];
                                            if (diagnostic_mode == 1)
                                            {
                                                Debug.WriteLine((row + i) + ", " + (col + j) + " is neighbour of depression , altidifference " + altidiff);
                                            }
                                            if (altidiff == minaltidiff)
                                            {   //if lowest higher nb = equally high as previous lowest higher nb
                                                lowestneighbourcounter++; numberoflowestneighbours++;
                                                if (numberoflowestneighbours == maxlowestnbs) { Debug.WriteLine(" WARNING: the setting for maximum number of lowest neighbours is " + numberoflowestneighbours + " lowest nbs" + maxlowestnbs); }
                                                rowlowestnb[lowestneighbourcounter] = (row + i);
                                                collowestnb[lowestneighbourcounter] = (col + j);
                                                // in this way, we can add all equally high lowest higher neighbours to the current depression (maximum = maxlowestnbs)
                                            } //end if higher neighbour, equal as before
                                            if (altidiff > minaltidiff || GlobalMethods.dtm[row + i, col + j] < depressionlevel[depressionnumber])  // this INCLUDES nodata cells bordering the lake!!
                                            {   //het hoogteverschil met deze buur telt alleen als minder hoog dan vorige buren OF lager dan meerniveau 
                                                minaltidiff = altidiff;
                                                if (GlobalMethods.dtm[row + i, col + j] < depressionlevel[depressionnumber] && already_lower_than_lakelevel == 1)
                                                {  // the new lowest neighbour is lower than lakelevel!! 
                                                    // We want to know all lowest nbs that are lower than lakelevel, so we do not zero the rowlowestnb 
                                                    // and colllowestnb arrays
                                                    lowestneighbourcounter++; numberoflowestneighbours++;
                                                    rowlowestnb[lowestneighbourcounter] = (row + i);
                                                    collowestnb[lowestneighbourcounter] = (col + j);
                                                    if (diagnostic_mode == 1) { Debug.WriteLine(" found another neighbour that is lower than lakelevel "); }
                                                    if (diagnostic_mode == 1) { Debug.WriteLine(" lower neighbour that is lower than lakelevel: " + rowlowestnb[lowestneighbourcounter] + ", " + collowestnb[lowestneighbourcounter] + " , " + GlobalMethods.dtm[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]]); }
                                                }
                                                else  // the new lowest neighbour is higher than lakelevel or is the first lowest nb that is lower than lakelevel
                                                {
                                                    if (GlobalMethods.dtm[row + i, col + j] < depressionlevel[depressionnumber]) { already_lower_than_lakelevel = 1; }
                                                    for (lowestneighbourcounter = 0; lowestneighbourcounter < maxlowestnbs; lowestneighbourcounter++)
                                                    {
                                                        rowlowestnb[lowestneighbourcounter] = -1; collowestnb[lowestneighbourcounter] = -1;
                                                    } //end:  for all lowestneighbours that we had before, the rowlowestnb and collowestnb arrays have been zeroed
                                                    lowestneighbourcounter = 0; numberoflowestneighbours = 1;
                                                    rowlowestnb[0] = (row + i); collowestnb[0] = (col + j);
                                                    if (diagnostic_mode == 1) { Debug.WriteLine(" higher neighbour that is lower than previous higher nbs: " + rowlowestnb[0] + ", " + collowestnb[0] + " , " + GlobalMethods.dtm[rowlowestnb[0], collowestnb[0]]); }
                                                }
                                            } //end if higher neighbour but lower than before
                                        } //end if nbismemberofdepression

                                    } // end if boundary
                                }
                            } // double end for circle with possibly extended radius around the sink, alpha - beta circle . We now know what is/are the lowest higher neighbours


                            if (diagnostic_mode == 1)
                            {
                                for (lowestneighbourcounter = 0; lowestneighbourcounter < numberoflowestneighbours; lowestneighbourcounter++)
                                {
                                    Debug.WriteLine(rowlowestnb[lowestneighbourcounter] + " " + collowestnb[lowestneighbourcounter] + " is one of the " + numberoflowestneighbours + " lowest neighbours of depression " + depressionnumber + ". Its altitude is " + GlobalMethods.dtm[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]]);
                                }
                            }

                            int outletfound = 0; int numberofoutlets = 0; int outletnumber = 0;

                            for (lowestneighbourcounter = 0; lowestneighbourcounter < numberoflowestneighbours; lowestneighbourcounter++)
                            {
                                if (rowlowestnb[lowestneighbourcounter] >= iupedge[depressionnumber]) { iupedge[depressionnumber] = rowlowestnb[lowestneighbourcounter] + 1; }
                                if (rowlowestnb[lowestneighbourcounter] <= iloedge[depressionnumber]) { iloedge[depressionnumber] = rowlowestnb[lowestneighbourcounter] - 1; }
                                if (collowestnb[lowestneighbourcounter] >= jupedge[depressionnumber]) { jupedge[depressionnumber] = collowestnb[lowestneighbourcounter] + 1; }
                                if (collowestnb[lowestneighbourcounter] <= jloedge[depressionnumber]) { jloedge[depressionnumber] = collowestnb[lowestneighbourcounter] - 1; }
                                if (diagnostic_mode == 1) { Debug.WriteLine(" minaltidiff = " + minaltidiff + "; depression " + depressionnumber + " row " + rowlowestnb[lowestneighbourcounter] + " col " + collowestnb[lowestneighbourcounter] + " alt " + GlobalMethods.dtm[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]] + " is lower than depressionlevel " + depressionlevel[depressionnumber]); }

                                if (minaltidiff <= 0.0)
                                { // if the cell is higher than depressionlevel
                                    // it can either be a cell of another depression

                                    if (diagnostic_mode == 1) { Debug.WriteLine(rowlowestnb[lowestneighbourcounter] + "," + collowestnb[lowestneighbourcounter] + " is member of depression " + depression[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]]); }
                                    if (depression[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]] != 0 && depression[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]] != depressionnumber)
                                    {  // then we have touched upon a depression that was previously analysed
                                        //GlobalMethods.status_map[rowlowestnb[lowestneighbourcounter],collowestnb[lowestneighbourcounter]] = 0;
                                        otherdepressionsize = 0;
                                        otherdepression = depression[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]];
                                        if (GlobalMethods.t > 1000000) { diagnostic_mode = 1; }
                                        totaldepressions--;
                                        totaldepressionvolume -= depressionvolume_m[otherdepression];
                                        for (int outletcounter = 0; outletcounter < 5; outletcounter++)
                                        {
                                            GlobalMethods.drainingoutlet_col[otherdepression, outletcounter] = -1;
                                            GlobalMethods.drainingoutlet_col[otherdepression, outletcounter] = -1;
                                        }
                                        depressionvolume_m[depressionnumber] += depressionvolume_m[otherdepression];
                                        //if (diagnostic_mode == 1) { Debug.WriteLine(" depressionvolume of depression " + depressionnumber + " was increased to " + depressionvolume[depressionnumber] + " with " + depressionvolume[otherdepression] + " of depression " + otherdepression); }
                                        depressionvolume_m[depressionnumber] += (depressionsize[depressionnumber] * (GlobalMethods.dtm[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]] - depressionlevel[depressionnumber]));  //
                                        //if (diagnostic_mode == 1) { Debug.WriteLine(" added " + depressionsize[depressionnumber] * (GlobalMethods.dtm[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]] - depressionlevel[depressionnumber]) + " to depressionvolume of depression " + depressionnumber + ". Dtm " + GlobalMethods.dtm[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]] + ", depressionlevel " + depressionlevel[depressionnumber]); }
                                        depressionvolume_m[otherdepression] = 0;
                                        depressionlevel[depressionnumber] = GlobalMethods.dtm[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]];
                                        //if (diagnostic_mode == 1) { Debug.WriteLine(" first: ilo = " + iloradius + " , iup = " + iupradius + " , jlo = " + jloradius + " , jup = " + jupradius + "  around sink " + row + " " + col); }
                                        //if (diagnostic_mode == 1) { Debug.WriteLine(" depression " + otherdepression + " : iloedge " + iloedge[otherdepression] + " , iupedge " + iupedge[otherdepression] + " , jloedge " + jloedge[otherdepression] + " , jupedge " + jupedge[otherdepression]); }
                                        if (jloedge[otherdepression] < (col - jloradius)) { jloradius = Math.Abs(jloedge[otherdepression] - col) + 1; jloedge[depressionnumber] = col - jloradius; }     // we enlarge the area around the sink of the current depression that is to be checked
                                        if (iloedge[otherdepression] < (row - iloradius)) { iloradius = Math.Abs(iloedge[otherdepression] - row) + 1; iloedge[depressionnumber] = row - iloradius; }     // it now includes the complete area of the depression that has been touched upon
                                        if (jupedge[otherdepression] > (col + jupradius)) { jupradius = Math.Abs(jupedge[otherdepression] - col) + 1; jupedge[depressionnumber] = col + jupradius; }
                                        if (iupedge[otherdepression] > (row + iupradius)) { iupradius = Math.Abs(iupedge[otherdepression] - row) + 1; iupedge[depressionnumber] = row + iupradius; }
                                        //if (diagnostic_mode == 1) { Debug.WriteLine(" now: ilo = " + iloradius + " , iup = " + iupradius + " , jlo = " + jloradius + " , jup = " + jupradius + " around sink " + row + "  " + col); }
                                        //if (diagnostic_mode == 1) { Debug.WriteLine(" begonnen met gebied om sink om cellen uit ander meer om te nummeren"); }
                                        for (alpha = (-1 * iloradius); alpha <= iupradius; alpha++)
                                        {       // move around in this square and change depressionnumber
                                            for (beta = (-1 * jloradius); beta <= jupradius; beta++)
                                            {
                                                if (((row + alpha) >= 0) && ((row + alpha) < GlobalMethods.nr) &&   // insofar that the circle is within the boundaries, excluding the centre cell itself
                                                    ((col + beta) >= 0) && ((col + beta) < GlobalMethods.nc) && !((alpha == 0) && (beta == 0)) && GlobalMethods.dtm[row + alpha, col + beta] != -9999)
                                                {
                                                    if (depression[row + alpha, col + beta] == otherdepression)
                                                    {
                                                        depression[row + alpha, col + beta] = depressionnumber;
                                                        otherdepressionsize++;
                                                        if (diagnostic_mode == 1) { Debug.WriteLine(" moved cell " + (row + alpha) + " , " + (col + beta) + "(" + GlobalMethods.dtm[row + alpha, col + beta] + ") from depression " + otherdepression + " to depression " + depressionnumber + " "); }
                                                    } // end if cell belonged to other depression
                                                }
                                            }
                                        } // double end for second alpha-beta square : around the sink to change depressionnumber of previously checked depression
                                        if (diagnostic_mode == 1) { Debug.WriteLine(" All " + otherdepressionsize + "  cells of depression " + otherdepression + "  were added to depression " + depressionnumber + "  (prvsly " + depressionsize[depressionnumber] + " ), level now " + depressionlevel[depressionnumber]); }
                                        depressionsize[depressionnumber] += otherdepressionsize;
                                        depressionsize[otherdepression] = 0;
                                        totaldepressionsize -= otherdepressionsize;
                                        if (diagnostic_mode == 1) { Debug.WriteLine("B totaldepressionsize " + totaldepressionsize + " , otherdepressionnsize " + otherdepression + "  = " + otherdepressionsize + "  depressionnumber " + depressionnumber + " "); }

                                        // it is theoretically possible that the outlet that connects the 'depression' and the 'otherdepression' , drains to a third side. In that case, the new, combined depression
                                        // should be declared ready after added 'otherdepression' to 'depression' ...

                                        for (alpha = (-1); alpha <= 1; alpha++)
                                        {
                                            for (beta = -1; beta <= 1; beta++)
                                            {
                                                if (((rowlowestnb[lowestneighbourcounter] + alpha) >= 0) && ((rowlowestnb[lowestneighbourcounter] + alpha) < GlobalMethods.nr) &&
                                                    ((collowestnb[lowestneighbourcounter] + beta) >= 0) && ((collowestnb[lowestneighbourcounter] + beta) < GlobalMethods.nc) &&
                                                    !((alpha == 0) && (beta == 0)) && GlobalMethods.dtm[rowlowestnb[lowestneighbourcounter] + alpha, collowestnb[lowestneighbourcounter] + beta] != -9999)
                                                {  // insofar that the circle is within the boundaries, excluding the centre cell itself
                                                    if (depression[rowlowestnb[lowestneighbourcounter] + alpha, collowestnb[lowestneighbourcounter] + beta] != depressionnumber)
                                                    {
                                                        if (GlobalMethods.dtm[rowlowestnb[lowestneighbourcounter] + alpha, collowestnb[lowestneighbourcounter] + beta] < depressionlevel[depressionnumber])
                                                        {
                                                            depressionready = 1;
                                                            GlobalMethods.drainingoutlet_col[depressionnumber, 0] = rowlowestnb[lowestneighbourcounter];
                                                            GlobalMethods.drainingoutlet_col[depressionnumber, 0] = collowestnb[lowestneighbourcounter];
                                                            GlobalMethods.status_map[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]] = 2;
                                                            if (diagnostic_mode == 1)
                                                            {
                                                                Debug.WriteLine(" depression " + depressionnumber + " drains to a third side and is ready ");
                                                                //displayonscreen(rowlowestnb[lowestneighbourcounter] + alpha, collowestnb[lowestneighbourcounter] + beta);
                                                            }
                                                        } // end if lower than outlet = depressionlevel
                                                    }  // end if depression != depression
                                                } // end if boundaries
                                            } // end for beta
                                        } // end for alpha
                                    }  // end if touched another depression with lowest nb

                                    //or it can be any other non-depression cell

                                    else
                                    {      // so we did not touch another depression with our lowest nb , but it was a higher or equally high nb so the depression is not yet ready
                                        if (diagnostic_mode == 1) { Debug.WriteLine(" this lowest neighbour: second option: no depression "); }
                                        depression[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]] = depressionnumber;
                                        if (diagnostic_mode == 1) { Debug.WriteLine(" depressionvolume of depression " + depressionnumber + " is " + depressionvolume_m[depressionnumber]); }
                                        depressionvolume_m[depressionnumber] += (depressionsize[depressionnumber] * (GlobalMethods.dtm[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]] - depressionlevel[depressionnumber]));  // add the amount of water added to the surface already part of depression
                                        if (diagnostic_mode == 1) { Debug.WriteLine(" added " + (depressionsize[depressionnumber] * (GlobalMethods.dtm[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]] - depressionlevel[depressionnumber])) + " to depressionvolume of depression " + depressionnumber + ". depressionsize: " + depressionsize[depressionnumber] + ", GlobalMethods.dtm %6.6f, depressionlevel " + GlobalMethods.dtm[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]], depressionlevel[depressionnumber]); }
                                        depressionsize[depressionnumber]++;
                                        depressionlevel[depressionnumber] = GlobalMethods.dtm[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]];
                                        if (diagnostic_mode == 1) { Debug.WriteLine(" added " + rowlowestnb[lowestneighbourcounter] + ", " + collowestnb[lowestneighbourcounter] + "  with level " + GlobalMethods.dtm[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]] + " to depression " + depressionnumber + " , size now " + depressionsize[depressionnumber] + " "); }
                                        if (diagnostic_mode == 1) { Debug.WriteLine(" lowest neighbour of depression " + depressionnumber + "  is " + rowlowestnb[lowestneighbourcounter] + " , " + collowestnb[lowestneighbourcounter] + " ," + depressionlevel[depressionnumber] + "  "); }
                                        if (diagnostic_mode == 1) { Debug.WriteLine(" depressionlevel for this depression is currently: " + depressionlevel[depressionnumber]); }  // laagste buur is nog niet betrokken bij een ander meer

                                        if (GlobalMethods.status_map[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]] == 3)
                                        { // then this depression drains out of the DEM...
                                            depressionready = 1;
                                            if (diagnostic_mode == 1) { Debug.WriteLine(" depression " + depressionnumber + " drains out of the DEM and is ready "); }
                                            depressiondrainsout = 1;
                                            GlobalMethods.drainingoutlet_col[depressionnumber, 0] = rowlowestnb[lowestneighbourcounter];
                                            GlobalMethods.drainingoutlet_col[depressionnumber, 0] = collowestnb[lowestneighbourcounter];
                                        }
                                        else
                                        { // depression does not drain out of DEM
                                            if ((rowlowestnb[lowestneighbourcounter] - row) == iupradius) { iupradius++; }     //in that case we will now change searchradius
                                            if ((row - rowlowestnb[lowestneighbourcounter]) == iloradius) { iloradius++; }
                                            if ((collowestnb[lowestneighbourcounter] - col) == jupradius) { jupradius++; }
                                            if ((col - collowestnb[lowestneighbourcounter]) == jloradius) { jloradius++; }
                                        } // end else
                                    } //end else
                                } // end if cell was higher than depressionlevel

                                else
                                {  // apparently it was lower than depressionlevel
                                    //Debug.WriteLine(" found lower neighbour " + rowlowestnb[lowestneighbourcounter] + "  " + collowestnb[lowestneighbourcounter] + "  alt " + GlobalMethods.dtm[rowlowestnb[lowestneighbourcounter], collowestnb[lowestneighbourcounter]] + "  for depression " + depressionnumber + "   level " + depressionlevel[depressionnumber] + " "); 
                                    outletfound = 1;
                                    // find for this cell, that should not be part of the depression, the necessarily present depression nb @ depressionlevel and call it a outlet
                                    for (alpha = (-1); alpha <= 1; alpha++)
                                    {
                                        for (beta = -1; beta <= 1; beta++)
                                        {
                                            if (((rowlowestnb[lowestneighbourcounter] + alpha) >= 0) && ((rowlowestnb[lowestneighbourcounter] + alpha) < GlobalMethods.nr) &&
                                                    ((collowestnb[lowestneighbourcounter] + beta) >= 0) && ((collowestnb[lowestneighbourcounter] + beta) < GlobalMethods.nc) &&
                                                    !((alpha == 0) && (beta == 0)) && GlobalMethods.dtm[rowlowestnb[lowestneighbourcounter] + alpha, collowestnb[lowestneighbourcounter] + beta] != -9999)
                                            {  // insofar that the circle is within the boundaries, excluding the centre cell itself
                                                if (depression[rowlowestnb[lowestneighbourcounter] + alpha, collowestnb[lowestneighbourcounter] + beta] == depressionnumber)
                                                {
                                                    if (GlobalMethods.dtm[rowlowestnb[lowestneighbourcounter] + alpha, collowestnb[lowestneighbourcounter] + beta] == depressionlevel[depressionnumber])
                                                    {
                                                        //Debug.WriteLine(" At lower cell  " + rowlowestnb[lowestneighbourcounter] + " " + collowestnb[lowestneighbourcounter] + ", looking back at " + (rowlowestnb[lowestneighbourcounter] + alpha) + " " + (collowestnb[lowestneighbourcounter] + beta));
                                                        if (GlobalMethods.drainingoutlet_col[depressionnumber, outletnumber] != -1) //ArT is outletnumber put to zero before here somewhere?
                                                        {
                                                            // Then we found an earlier outlet (somewhere), and we have to define an extra one here, IF that earlier one was not this one.      
                                                            bool this_is_an_earlier_outlet = false;
                                                            for (int outletcounter = 0; outletcounter < 5; outletcounter++)
                                                            {
                                                                if (rowlowestnb[lowestneighbourcounter] + alpha == GlobalMethods.drainingoutlet_col[depressionnumber, outletcounter] && collowestnb[lowestneighbourcounter] + beta == GlobalMethods.drainingoutlet_col[depressionnumber, outletcounter])
                                                                {
                                                                    this_is_an_earlier_outlet = true;
                                                                    //Debug.WriteLine(" A prior outlet definition for depression " + depressionnumber + " exists at " + GlobalMethods.drainingoutlet_col[depressionnumber, outletcounter] + " " + GlobalMethods.drainingoutlet_col[depressionnumber, outletcounter] + " (" + (lowestneighbourcounter + 1) + "/" + numberoflowestneighbours + ")");
                                                                }
                                                            }
                                                            if (this_is_an_earlier_outlet == false)
                                                            {
                                                                //Debug.WriteLine(" A new outlet was found for depression " + depressionnumber + " at " + (rowlowestnb[lowestneighbourcounter] + alpha) + " " + (collowestnb[lowestneighbourcounter] + beta) + " (" + (lowestneighbourcounter + 1) + "/" + numberoflowestneighbours + ")");
                                                                outletnumber++;
                                                                numberofoutlets++;
                                                                if (outletnumber > 4)
                                                                {
                                                                    //displayonscreen(GlobalMethods.drainingoutlet_col[depressionnumber, 1], GlobalMethods.drainingoutlet_col[depressionnumber, 1]);
                                                                    if (diagnostic_mode == 1) { Debug.WriteLine(" Warning: LORICA seeks to define more than five outlets for depression " + depressionnumber + ". This request is denied"); }
                                                                    outletnumber--; numberofoutlets--;
                                                                }
                                                                else
                                                                {
                                                                    GlobalMethods.drainingoutlet_col[depressionnumber, outletnumber] = rowlowestnb[lowestneighbourcounter] + alpha;
                                                                    GlobalMethods.drainingoutlet_col[depressionnumber, outletnumber] = collowestnb[lowestneighbourcounter] + beta;
                                                                    GlobalMethods.status_map[rowlowestnb[lowestneighbourcounter] + alpha, collowestnb[lowestneighbourcounter] + beta] = 2;
                                                                    //Debug.WriteLine(" made " + (rowlowestnb[lowestneighbourcounter] + alpha) + "  " + (collowestnb[lowestneighbourcounter] + beta) + "  the GlobalMethods.nr " + outletnumber + " outlet for depression " + depressionnumber + " ");
                                                                }
                                                            }

                                                        }
                                                        else  // then this is the first outlet, and we will define it as such
                                                        {
                                                            //Debug.WriteLine(" First definition of outlet for depression " + depressionnumber + " at " + (rowlowestnb[lowestneighbourcounter] + alpha) + " " + (collowestnb[lowestneighbourcounter] + beta) + " (" + (lowestneighbourcounter+1) + "/" + numberoflowestneighbours + ")");
                                                            if (outletnumber > 4)
                                                            {
                                                                //displayonscreen(GlobalMethods.drainingoutlet_col[depressionnumber, 1], GlobalMethods.drainingoutlet_col[depressionnumber, 1]);
                                                                if (diagnostic_mode == 1) { Debug.WriteLine(" Warning: LORICA seeks to define more than five outlets for depression " + depressionnumber + ". This request is denied"); }
                                                                outletnumber--; numberofoutlets--;
                                                            }
                                                            else
                                                            {
                                                                GlobalMethods.drainingoutlet_col[depressionnumber, outletnumber] = rowlowestnb[lowestneighbourcounter] + alpha;
                                                                GlobalMethods.drainingoutlet_col[depressionnumber, outletnumber] = collowestnb[lowestneighbourcounter] + beta;
                                                                GlobalMethods.status_map[rowlowestnb[lowestneighbourcounter] + alpha, collowestnb[lowestneighbourcounter] + beta] = 2;
                                                                //Debug.WriteLine(" made " + (rowlowestnb[lowestneighbourcounter] + alpha) + "  " + (collowestnb[lowestneighbourcounter] + beta) + "  the GlobalMethods.nr " + outletnumber + " outlet for depression " + depressionnumber + " ");
                                                            }
                                                        }

                                                    } // end if GlobalMethods.dtm = depressionlevel
                                                } // end if depression = depressionnumber
                                            } // end if bnd
                                        } // end for beta
                                    } // end for alpha
                                } // end else
                            } //end for all lowestneighbours
                            // if outlet(s) defined ==> depressionready
                            if (outletfound == 1) { depressionready = 1; }
                        } // end while depressionready != 1

                        //if (depressionnumber == 180) {diagnostic_mode = 1;}
                        if (diagnostic_mode == 1)
                        {
                            Debug.WriteLine(" depression " + depressionnumber + "  (" + depressionlevel[depressionnumber] + "  m) is ready and contains " + depressionsize[depressionnumber] + "  cells. Volume = " + depressionvolume_m[depressionnumber] + " ");
                            if (depressiondrainsout == 1) { Debug.WriteLine(" depression " + depressionnumber + "  (" + depressionlevel[depressionnumber] + "  m) drains outside the DEM "); }
                            minimaps(row, col);
                        }
                        totaldepressionsize += depressionsize[depressionnumber];
                        totaldepressionvolume += depressionvolume_m[depressionnumber];
                        if (maxsize < depressionsize[depressionnumber]) { maxsize = depressionsize[depressionnumber]; largestdepression = depressionnumber; }
                        if (maxdepressionnumber < depressionnumber) { maxdepressionnumber = depressionnumber; }
                        if (diagnostic_mode == 1)
                        {
                            Debug.WriteLine(" Defined depression " + depressionnumber + ". Now size " + depressionsize[depressionnumber] + " and volume " + depressionvolume_m[depressionnumber]);
                        }
                    } // end if sink
                    //Debug.WriteLine("now at row " + row + " and col " + col);
                }  // end for  col
            } // end for   row

            Debug.WriteLine("\n\n--depressions overview--");
            if (totaldepressions != 0)
            {
                Debug.WriteLine("found " + totaldepressions + "  depressions containing " + totaldepressionsize + "  cells, with a volume of " + totaldepressionvolume);
                Debug.WriteLine(" " + totaldepressions + " depressions with a volume of " + totaldepressionvolume);
                //Debug.WriteLine("depression " + largestdepression + "  is largest by area: " + maxsize + " cells " + depressionlevel[largestdepression] + " m " + depressionvolume[largestdepression] + "m3");
                //if (depressionvolume[largestdepression] < 0) { Debug.WriteLine("negative depressionvolume found"); }
            }
            else
            {
                Debug.WriteLine(" no depressions found ");

            }
            //Debug.WriteLine("\n");
            //Debug.WriteLine(" drains at: " + + " , " + + " , volume: %6.9f " ,GlobalMethods.drainingoutlet_col[largestdepression],GlobalMethods.drainingoutlet_col[largestdepression],depressionvolume[largestdepression]); 
            //out_integer("lakes.asc",depression);
            //out_integer("status.asc", GlobalMethods.status_map);
            // 
            // * */
        }

        void define_fillheight_new()  //calculates where depressions must be filled how high
        {
            // when completely filling a depression, we need to know - for each constituent cell and even for its neighbours - the altitude we
            // can fill it to. Since the depression must still drain towards the outlet, we add a very small value
            // to membercells so they drain towards the outlet.
            // we cannot simply use distance_to_outlet for each member cell, since depressions can round corners....

            this.InfoStatusPanel.Text = "def fillheight has been entered";
            //Debug.WriteLine("defining fillheight\n");
            int notyetdone, done, depressiontt;

            once_dtm_fill = 0;

            for (row = 0; row < GlobalMethods.nr; row++)
            {
                for (col = 0; col < GlobalMethods.nc; col++)
                {
                    GlobalMethods.dtmfill_A[row, col] = -1;
                } //for
            } //for


            depressiontt = 0;
            for (depressiontt = 1; depressiontt <= maxdepressionnumber; depressiontt++)
            {  // for all possible depressions

                if (depressionsize[depressiontt] > 0)
                {      // if they exist, so have not been intermediate depressions
                    GlobalMethods.dtmfill_A[GlobalMethods.drainingoutlet_col[depressiontt, 0], GlobalMethods.drainingoutlet_col[depressiontt, 0]] = depressionlevel[depressiontt];
                    if (GlobalMethods.drainingoutlet_col[depressiontt, 1] != -1) { GlobalMethods.dtmfill_A[GlobalMethods.drainingoutlet_col[depressiontt, 1], GlobalMethods.drainingoutlet_col[depressiontt, 1]] = depressionlevel[depressiontt]; }
                    if (GlobalMethods.drainingoutlet_col[depressiontt, 2] != -1) { GlobalMethods.dtmfill_A[GlobalMethods.drainingoutlet_col[depressiontt, 2], GlobalMethods.drainingoutlet_col[depressiontt, 2]] = depressionlevel[depressiontt]; }
                    if (GlobalMethods.drainingoutlet_col[depressiontt, 3] != -1) { GlobalMethods.dtmfill_A[GlobalMethods.drainingoutlet_col[depressiontt, 3], GlobalMethods.drainingoutlet_col[depressiontt, 3]] = depressionlevel[depressiontt]; }
                    if (GlobalMethods.drainingoutlet_col[depressiontt, 4] != -1) { GlobalMethods.dtmfill_A[GlobalMethods.drainingoutlet_col[depressiontt, 4], GlobalMethods.drainingoutlet_col[depressiontt, 4]] = depressionlevel[depressiontt]; }

                    notyetdone = 1; done = 0;
                    while (notyetdone > 0)
                    {
                        notyetdone = 0;
                        //if (diagnostic_mode == 1) { Debug.WriteLine("depressioncells depression " + depressiontt + " size " + depressionsize[depressiontt] + " " +GlobalMethods.drainingoutlet_col[depressiontt]+ " " + GlobalMethods.drainingoutlet_col[depressiontt]); }
                        //diagnostic_mode = 0;
                        for (row = iloedge[depressiontt]; row <= iupedge[depressiontt]; row++)
                        {
                            for (col = jloedge[depressiontt]; col <= jupedge[depressiontt]; col++)
                            {
                                if (((row) >= 0) && ((row) < GlobalMethods.nr) && ((col) >= 0) && ((col) < GlobalMethods.nc) && GlobalMethods.dtm[row, col] != -9999)
                                {  //bnd
                                    if (GlobalMethods.t > 1000000) { diagnostic_mode = 1; } else { diagnostic_mode = 0; }
                                    if (diagnostic_mode == 1) { Debug.WriteLine("GlobalMethods.dtmfill_A of " + row + " " + col + " = " + GlobalMethods.dtmfill_A[row, col] + ", checking on behalf of depression " + depressiontt); }
                                    if (GlobalMethods.dtmfill_A[row, col] == -1 && depression[row, col] == depressiontt)
                                    {  // if this is a cell of the current depression that has not yet got a dtmfill
                                        // and remember that this is the only place where this cell could have gotten that DTMfill
                                        notyetdone++;       // then we are not yet ready
                                        if (diagnostic_mode == 1) { Debug.WriteLine(row + " " + col + ": notyetdone =  " + notyetdone); }
                                        //if (diagnostic_mode == 1) { displayonscreen(row, col); }
                                        for (i = -1; i <= 1; i++)
                                        {      // go and see if it has a depression-nb that does have a dtmfill and that is lower than the previous possible dep-nb's dtmfill
                                            for (j = -1; j <= 1; j++)
                                            {
                                                if (((row + i) >= 0) && ((row + i) < GlobalMethods.nr) && ((col + j) >= 0) && ((col + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)) && GlobalMethods.dtm[row + i, col + j] != -9999)
                                                {  //bnd
                                                    if (depression[row + i, col + j] == depressiontt && GlobalMethods.dtmfill_A[row + i, col + j] > 0.0)
                                                    { // if it IS a depression-nb and DOES have a dtmfill
                                                        if ((i == 0 || j == 0) && (GlobalMethods.dtmfill_A[row, col] == -1 || (GlobalMethods.dtmfill_A[row, col] > GlobalMethods.dtmfill_A[row + i, col + j] + epsilon * GlobalMethods.dx)))
                                                        {
                                                            GlobalMethods.dtmfill_A[row, col] = (GlobalMethods.dtmfill_A[row + i, col + j] + epsilon * GlobalMethods.dx);
                                                        } // end if
                                                        if ((i != 0 && j != 0) && (GlobalMethods.dtmfill_A[row, col] == -1 || (GlobalMethods.dtmfill_A[row, col] > GlobalMethods.dtmfill_A[row + i, col + j] + epsilon * GlobalMethods.dx * Math.Sqrt(2))))
                                                        {
                                                            GlobalMethods.dtmfill_A[row, col] = (GlobalMethods.dtmfill_A[row + i, col + j] + epsilon * GlobalMethods.dx * Math.Sqrt(2));
                                                        } // end if
                                                    } // end if
                                                }//end if within boundaries
                                            }//end for circle (first part)
                                        }//end for circle (second part)
                                        if (GlobalMethods.dtmfill_A[row, col] > 0) { notyetdone--; done++; }
                                        //if (diagnostic_mode == 1) { Debug.WriteLine(" notyetdone =  " + notyetdone + " done = " + done); }
                                    } // end if depression
                                    //if (diagnostic_mode == 1) { Debug.WriteLine("-> " + notyetdone); }
                                } // end if bnd
                            } // end for
                        } // end for
                    } // end while
                } // end if they exist
            } //end for all possible depressions

            Debug.WriteLine("\n--dtmfill determination finished--");

        }

        void cleardelta(int iloradius, int iupradius, int jloradius, int jupradius, int clear_row, int clear_col)   //clears a delta
        {
            int epsilon, eta;
            if (diagnostic_mode == 1) { Debug.WriteLine(" clearing delta lake " + Math.Abs(depression[clear_row, clear_col]) + " around " + clear_row + " " + clear_col); }
            for (epsilon = -(iloradius + 3); epsilon <= iupradius + 3; epsilon++)
            {
                for (eta = -(jloradius + 3); eta <= jupradius + 3; eta++)
                {
                    if (((clear_row + epsilon) >= 0) && ((clear_row + epsilon) < GlobalMethods.nr) && ((clear_col + eta) >= 0) && ((clear_col + eta) < GlobalMethods.nc))
                    { // boundaries
                        if (depression[clear_row + epsilon, clear_col + eta] < 0)
                        {
                            depression[clear_row + epsilon, clear_col + eta] = Math.Abs(depression[clear_row, clear_col]);
                            if (diagnostic_mode == 1) { Debug.WriteLine(" membership of delta has been cancelled for " + (clear_row + epsilon) + " " + (clear_col + eta)); }
                        }  // end if depression < 0
                    }   // end if boundary
                }   // end for eta
            } // end for epsilon
            //Debug.WriteLine("cleared delta\n");

        }

        void update_depression(int number)   //updates depressions when the erosion/deposition process has reached them to include cells that have been eroded to below lakelevel
        {
            /*  a.	First estimate of required sediment is fillheight – GlobalMethods.dtm for all lake cells. Test whether GlobalMethods.dz_ero_m and GlobalMethods.dz_sed_m for these cells are zero (they should be).
                b.	Starting from every lake cell, look for cells around it that are not part of the lake and had GlobalMethods.dtm above lakelevel.
                    •	For such a cell, look around it for all lake-neighbours and determine the one that would yield the lowest fillheight. 
                    •	Assign that fillheight to the cell, see  if its current altitude (corrected for already calculated ero and sed) is lower than fillheight.
                    •	Add the difference to the sediment needed to fill the (now larger) lake.
                    •	Add the cell to the lake and update the size and volume  of the lake
                c.	Continue until no more cells around the lake are lower 
                This would potentially clash with lakes that have more than two outlets. The third and later outlets would be 
                seen as lake cells and their lower neighbours on the non-lake-side would be added to the lake, leading to errors. */
            int urow = 0, ucol = 0, size = 0;
            int depressionnumber = number;
            //if (depressionnumber > 1) diagnostic_mode = 1;
            if (diagnostic_mode == 1) { Debug.WriteLine(" now updating depression " + depressionnumber); }
            depressionsum_water_m = 0;
            depressionsum_sediment_m = 0;
            depressionsum_texture_kg[0] = 0; depressionsum_texture_kg[1] = 0; depressionsum_texture_kg[2] = 0; depressionsum_texture_kg[3] = 0; depressionsum_texture_kg[4] = 0; depressionsum_YOM_kg = 0; depressionsum_OOM_kg = 0;
            needed_to_fill_depression_m = 0;

            for (urow = iloedge[depressionnumber]; urow <= iupedge[depressionnumber]; urow++)
            {
                for (ucol = jloedge[depressionnumber]; ucol <= jupedge[depressionnumber]; ucol++)
                {
                    if (((urow) >= 0) && ((urow) < GlobalMethods.nr) && ((ucol) >= 0) && ((ucol) < GlobalMethods.nc) && GlobalMethods.dtm[urow, ucol] != -9999)
                    {  //bnd
                        if (depression[urow, ucol] == depressionnumber)
                        {
                            if (only_waterflow_checkbox.Checked == false)
                            {
                                if (!(GlobalMethods.dz_ero_m[urow, ucol] == 0 && GlobalMethods.dz_sed_m[urow, ucol] == 0))
                                {
                                    Debug.WriteLine(" error in depression " + depressionnumber);
                                    // this should not happen: erosion or sedimentation into lake cells is not allowed - only the provision of those cells with sediment in transport.
                                    minimaps(urow, ucol);
                                }
                            }
                            // if they are member of the lake : add the diff between fillheight and GlobalMethods.dtm to the volume that must be filled, and add the sedintrans to what is available for filling
                            depressionsum_water_m += GlobalMethods.waterflow_m3[urow, ucol] / GlobalMethods.dx / GlobalMethods.dx;
                            if (only_waterflow_checkbox.Checked == false)
                            {
                                for (size = 0; size < GlobalMethods.n_texture_classes; size++)
                                {
                                    depressionsum_texture_kg[size] += GlobalMethods.sediment_in_transport_kg[urow, ucol, size];
                                }
                                depressionsum_OOM_kg += GlobalMethods.old_SOM_in_transport_kg[urow, ucol];
                                depressionsum_YOM_kg += GlobalMethods.young_SOM_in_transport_kg[urow, ucol];
                            }
                            needed_to_fill_depression_m += GlobalMethods.dtmfill_A[urow, ucol] - GlobalMethods.dtm[urow, ucol];
                            //if (diagnostic_mode == 1) { Debug.WriteLine(" dep_cell " + row + " " + col + " fillheight " + GlobalMethods.dtmfill_A[urow, ucol] + " GlobalMethods.dtm " + GlobalMethods.dtm[urow, ucol] + " needed now " + needed_to_fill_depression  ); }
                        }
                    }
                }
            }

            //now that we know how much kgs of every material are available for lake filling, we can calculate how much thickness [m3/m2 = m] that means.
            depressionsum_sediment_m = calc_thickness_from_mass(depressionsum_texture_kg, depressionsum_YOM_kg, depressionsum_OOM_kg);


            int updating_lake = 1;
            while (updating_lake == 1)
            {
                // while there are potential cells to be added to the lake
                updating_lake = 0;
                for (urow = iloedge[depressionnumber] - 1; urow <= iupedge[depressionnumber] + 1; urow++)
                {
                    for (ucol = jloedge[depressionnumber] - 1; ucol <= jupedge[depressionnumber] + 1; ucol++)
                    {
                        if (((urow) >= 0) && ((urow) < GlobalMethods.nr) && ((ucol) >= 0) && ((ucol) < GlobalMethods.nc) && GlobalMethods.dtm[urow, ucol] != -9999)
                        {  //bnd
                            if (depression[urow, ucol] != depressionnumber && GlobalMethods.dtm[urow, ucol] > depressionlevel[depressionnumber]) // the second part of the condition should ensure that no cells on the downstream side of outlets are added. 
                            {
                                double lowest_dtm_fill = 9999;
                                // we have now found a cell that should potentiall added to the lake - however, it should have a lake-neighbour for that to be really true. We look around to check this.
                                for (alpha = -1; alpha <= 1; alpha++)
                                {
                                    for (beta = -1; beta <= 1; beta++)
                                    {
                                        if ((urow + alpha) >= 0 && (urow + alpha) < GlobalMethods.nr && (ucol + beta) >= 0 && (ucol + beta) < GlobalMethods.nc && !((alpha) == 0 && (beta) == 0) && GlobalMethods.dtm[urow + alpha, ucol + beta] != -9999)
                                        {
                                            if ((urow != urow + alpha) && (ucol != ucol + beta)) { GlobalMethods.d_x = GlobalMethods.dx * Math.Sqrt(2); } else { GlobalMethods.d_x = GlobalMethods.dx; }
                                            if (only_waterflow_checkbox.Checked == false)
                                            {
                                                if (depression[urow + alpha, ucol + beta] == depressionnumber && (GlobalMethods.dtm[urow, ucol] + GlobalMethods.dz_ero_m[urow, ucol] + GlobalMethods.dz_sed_m[urow, ucol]) < (GlobalMethods.dtmfill_A[urow + alpha, ucol + beta] + epsilon * GlobalMethods.d_x))
                                                // if this cell has a lake neighbour and has a current altitude lower than the fillheight resulting from that neighbour, we store that fillheight until we find the lowest fillheight resulting from any lake-nbs of this cell
                                                {
                                                    if ((GlobalMethods.dtmfill_A[urow + alpha, ucol + beta] + epsilon * GlobalMethods.d_x) < lowest_dtm_fill)
                                                    {
                                                        lowest_dtm_fill = (GlobalMethods.dtmfill_A[urow + alpha, ucol + beta] + epsilon * GlobalMethods.d_x);

                                                    }
                                                }
                                            }
                                            else
                                            {
                                                if (depression[urow + alpha, ucol + beta] == depressionnumber && (GlobalMethods.dtm[urow, ucol]) < (GlobalMethods.dtmfill_A[urow + alpha, ucol + beta] + epsilon * GlobalMethods.d_x))
                                                // if this cell has a lake neighbour and has a current altitude lower than the fillheight resulting from that neighbour, we store that fillheight until we find the lowest fillheight resulting from any lake-nbs of this cell
                                                {
                                                    if ((GlobalMethods.dtmfill_A[urow + alpha, ucol + beta] + epsilon * GlobalMethods.d_x) < lowest_dtm_fill)
                                                    {
                                                        lowest_dtm_fill = (GlobalMethods.dtmfill_A[urow + alpha, ucol + beta] + epsilon * GlobalMethods.d_x);

                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                // after checking, if this cell had any lake-nb and was low enough itself, it is added to the lake and the search window is enlarged.
                                if (lowest_dtm_fill < 9999)
                                {
                                    depression[urow, ucol] = depressionnumber;
                                    depressionsize[depressionnumber]++;
                                    GlobalMethods.dtmfill_A[urow, ucol] = lowest_dtm_fill;
                                    if (only_waterflow_checkbox.Checked == false)
                                    {
                                        needed_to_fill_depression_m += GlobalMethods.dtmfill_A[urow, ucol] - (GlobalMethods.dtm[urow, ucol] + GlobalMethods.dz_ero_m[urow, ucol] + GlobalMethods.dz_sed_m[urow, ucol]);
                                        for (size = 0; size < GlobalMethods.n_texture_classes; size++)
                                        {
                                            GlobalMethods.sediment_in_transport_kg[urow, ucol, size] = 0;
                                        }
                                    }
                                    else { needed_to_fill_depression_m += GlobalMethods.dtmfill_A[urow, ucol] - (GlobalMethods.dtm[urow, ucol]); }
                                    updating_lake = 1;
                                    if (urow == iloedge[depressionnumber] - 1) { iloedge[depressionnumber]--; }
                                    if (urow == iupedge[depressionnumber] + 1) { iupedge[depressionnumber]++; }
                                    if (ucol == jloedge[depressionnumber] - 1) { jloedge[depressionnumber]--; }
                                    if (ucol == jupedge[depressionnumber] + 1) { jupedge[depressionnumber]++; }

                                    // setting sed in trans to 0 is required to avoid double-counting when building delta's later
                                    // no addition to depressionsum_sed or _water is required because flows from this cell have already been considered earlier, and any GlobalMethods.dz_ero_m or GlobalMethods.dz_sed_m have been considered, and sed_in_trans has arrived at the lake to be counted
                                    //if (diagnostic_mode == 1) { Debug.WriteLine(" non_dep_cell " + urow + " " + ucol + " fillheight " + GlobalMethods.dtmfill_A[urow, ucol] + " new GlobalMethods.dtm " + (GlobalMethods.dtm[urow, ucol] + GlobalMethods.dz_ero_m[urow, ucol] + GlobalMethods.dz_sed_m[urow, ucol]) + " needed now " + needed_to_fill_depression); }
                                }
                            }
                        }
                    }
                }  // end for

            } // end while
            if (diagnostic_mode == 1) { Debug.WriteLine(" Updated depression " + depressionnumber + ". Now size " + depressionsize[depressionnumber] + ", sum_water " + depressionsum_water_m + " sum_sed " + depressionsum_sediment_m + " Needed: " + needed_to_fill_depression_m); }

        }

        void fill_depression(int number)  // completely fills updated depressions (because enough sediment was available)
        {
            int this_depression = number;
            int fillrow = 0, size, fillcol = 0;
            //Debug.WriteLine(" now filling depression " + this_depression);
            sediment_filled_m += needed_to_fill_depression_m;


            for (fillrow = iloedge[this_depression]; fillrow <= iupedge[this_depression]; fillrow++)
            {
                for (fillcol = jloedge[this_depression]; fillcol <= jupedge[this_depression]; fillcol++)
                {
                    if (((fillrow) >= 0) && ((fillrow) < GlobalMethods.nr) && ((fillcol) >= 0) && ((fillcol) < GlobalMethods.nc) && GlobalMethods.dtm[fillrow, fillcol] != -9999)
                    {  //bnd
                        if (depression[fillrow, fillcol] == this_depression)
                        {
                            //Debug.WriteLine(" adding "  +(GlobalMethods.dtmfill_A[fillrow, fillcol] - GlobalMethods.dz_ero_m[fillrow, fillcol] - GlobalMethods.dz_sed_m[fillrow, fillcol] - GlobalMethods.dtm[fillrow, fillcol])+" sed to " + fillrow + " "+ fillcol ); 
                            //sediment_filled2 += GlobalMethods.dtmfill_A[fillrow, fillcol] - GlobalMethods.dz_ero_m[fillrow, fillcol] - GlobalMethods.dz_sed_m[fillrow, fillcol] - GlobalMethods.dtm[fillrow, fillcol];   //sediment_filled2 SHOULD BE sediment_filled at the end
                            GlobalMethods.lake_sed_m[fillrow, fillcol] += GlobalMethods.dtmfill_A[fillrow, fillcol] - GlobalMethods.dz_ero_m[fillrow, fillcol] - GlobalMethods.dz_sed_m[fillrow, fillcol] - GlobalMethods.dtm[fillrow, fillcol];
                            GlobalMethods.dtm[fillrow, fillcol] = GlobalMethods.dtmfill_A[fillrow, fillcol] - GlobalMethods.dz_ero_m[fillrow, fillcol] - GlobalMethods.dz_sed_m[fillrow, fillcol];
                            if (GlobalMethods.dtm[fillrow, fillcol] == -1) { Debug.WriteLine("C cell " + (fillrow) + " " + (fillcol) + " has an altitude of -1 now"); minimaps(fillrow, fillcol); }
                        }
                    }
                }
            }
            int outletcounter = 0;
            while (GlobalMethods.drainingoutlet_col[this_depression, outletcounter] != -1)
            {
                outletcounter++;
                if (outletcounter == 5) { break; }
            }
            for (i = 0; i < outletcounter; i++)
            {
                GlobalMethods.waterflow_m3[GlobalMethods.drainingoutlet_col[this_depression, i], GlobalMethods.drainingoutlet_col[this_depression, i]] += GlobalMethods.dx * GlobalMethods.dx * (depressionvolume_m[this_depression]) / outletcounter;
                //we need a rule to determine whether any sediment size stays behind preferentially if we fill a depression. For the moment: no preference
                for (size = 0; size < GlobalMethods.n_texture_classes; size++)
                {
                    //GlobalMethods.sediment_in_transport_kg[GlobalMethods.drainingoutlet_col[this_depression, i], GlobalMethods.drainingoutlet_col[this_depression, i]] += (depressionsum_sediment_m - needed_to_fill_depression_m) / outletcounter * ;
                    //we need to calculate how much of which fraction there is in transport, and then deposit by weight-ratio (or anything else).                
                }
            }
            //Debug.WriteLine(" Filled depression " + this_depression);
        }

        void leave_depression_alone(int number)  // only updates counters and sentinels
        {
            /*int this_depression = number;
            int leaverow = 0, leavecol = 0;
            //Debug.WriteLine(" leaving depression " + this_depression + " alone");

            for (leaverow = iloedge[this_depression];
 leaverow <= iupedge[this_depression]; leaverow++)
            {
                for (leavecol = jloedge[this_depression]; leavecol <= jupedge[this_depression]; leavecol++)
                {
                    if (((leaverow) >= 0) && ((leaverow) < GlobalMethods.nr) && ((leavecol) >= 0) && ((leavecol) < GlobalMethods.nc) && GlobalMethods.dtm[leaverow, leavecol] != -9999)
                    {  //bnd
                        if (depression[leaverow, leavecol] == this_depression)
                        {
                            considered[leaverow, leavecol] = 1;
                        }
                    }
                }
            }
            int outletcounter = 0;
            while (GlobalMethods.drainingoutlet_col[this_depression, outletcounter] != -1)
            {
                outletcounter++;
                if (outletcounter == 5) { break; }
            }
            for (i = 0; i < outletcounter; i++)
            {
                considered[GlobalMethods.drainingoutlet_col[this_depression, i], GlobalMethods.drainingoutlet_col[this_depression, i]] = 0;
            }
            //Debug.WriteLine(" Left depression " + this_depression + " alone"); */
        }

        void delta_depression(int number)  // builds deltas in an updated depression (because not enough sed)
        {
            /*When there is not enough sediment to fill the lake, fill the lake from each of its initial side-cells that have a non-zero sediment in transport, excluding outlet cells (their sediment in transport gets moved outside – they will never have to be raised higher than lakelevel which they already have  - they will remain having the lakenumber even when the whole lake would have been filled). 

                c. For each of these cells, while there is sediment in transport /GlobalMethods.dx left:
                    •	Find a potentially lower oblique cell (in the lake, but there are  no others)
                    •	Find the first higher oblique cell relative to that cell
                    •	Raise the oblique deepest cell to that level
                    •	Reduce the remaining amount of sediment in transport with the raised amount
                    •	Test whether the cell has now been raised above fillheight, in which case: 
                            it must be lowered to its fillheight, 
                            remaining amount of sed in trans must be increased again
                            lakecell must be removed from the lake (potentially fragmenting the original lake) */

            int active_depression = number, size;
            if (diagnostic_mode == 1) { Debug.WriteLine(" building deltas in depression " + number + " sed needed " + needed_to_fill_depression_m + " sed available " + depressionsum_sediment_m); }


            //if (number == 30001 && GlobalMethods.t == 1) { diagnostic_mode = 1; } else { diagnostic_mode = 0; }

            for (startrow = iloedge[active_depression]; startrow <= iupedge[active_depression]; startrow++)
            {
                for (startcol = jloedge[active_depression]; startcol <= jupedge[active_depression]; startcol++)
                {
                    if (((startrow) >= 0) && ((startrow) < GlobalMethods.nr) && ((startcol) >= 0) && ((startcol) < GlobalMethods.nc) && GlobalMethods.dtm[startrow, startcol] != -9999)
                    {  //bnd
                        int sediment_present = 0;
                        for (size = 0; size < GlobalMethods.n_texture_classes; size++)
                        {
                            if (GlobalMethods.sediment_in_transport_kg[startrow, startcol, size] > 0)
                            { sediment_present = 1; }
                        }
                        if (depression[startrow, startcol] == active_depression && sediment_present == 1)
                        // so, for all cells in the depression that have a non-zero sediment in transport (and EXCLUDING the outlets) 
                        {
                            if (!(GlobalMethods.drainingoutlet_col[active_depression, 0] == startrow && GlobalMethods.drainingoutlet_col[active_depression, 0] == startcol) &&
                                !(GlobalMethods.drainingoutlet_col[active_depression, 1] == startrow && GlobalMethods.drainingoutlet_col[active_depression, 1] == startcol) &&
                                !(GlobalMethods.drainingoutlet_col[active_depression, 2] == startrow && GlobalMethods.drainingoutlet_col[active_depression, 2] == startcol) &&
                                !(GlobalMethods.drainingoutlet_col[active_depression, 3] == startrow && GlobalMethods.drainingoutlet_col[active_depression, 3] == startcol) &&
                                !(GlobalMethods.drainingoutlet_col[active_depression, 4] == startrow && GlobalMethods.drainingoutlet_col[active_depression, 4] == startcol))
                            {
                                deltasize = 0;
                                dhobliquemax1 = 0;
                                iloradius3 = 1; iupradius3 = 1; jloradius3 = 1; jupradius3 = 1;
                                iloradius2 = 1; iupradius2 = 1; jloradius2 = 1; jupradius2 = 1;
                                rowlowestobnb = startrow; collowestobnb = startcol;
                                if (diagnostic_mode == 1) { Debug.WriteLine(" building a delta from " + startrow + " " + startcol + ", GlobalMethods.dtm " + GlobalMethods.dtm[startrow, startcol]); minimaps(startrow, startcol); }
                                double[] local_s_i_t_kg = new double[5] { 0, 0, 0, 0, 0 };
                                for (size = 0; size < GlobalMethods.n_texture_classes; size++)
                                {
                                    local_s_i_t_kg[size] = GlobalMethods.sediment_in_transport_kg[startrow, startcol, size] / GlobalMethods.dx;
                                    GlobalMethods.sediment_in_transport_kg[startrow, startcol, size] = 0;
                                }
                                available_for_delta_m += calc_thickness_from_mass(local_s_i_t_kg, 0, 0);
                                sediment_delta_m += available_for_delta_m;
                                // we will completely use all of this sed in trans now, so let's set it to zero already
                                // and add it to the overall counter of the volume used in delta's

                                while (available_for_delta_m > 0)
                                {

                                    find_lowest_oblique_neighbour(active_depression);  // to start in the right location with building the delta
                                    // this may be lower in the lake than at the side-cell, at location rowlowestobnb,collowestobnb
                                    depression[rowlowestobnb, collowestobnb] = -active_depression;   // we have found the start of this delta
                                    deltasize = 1;    // therefore we give delta a value of 1
                                    dhobliquemax2 = 99999.99;
                                    II = 0; JJ = 0;
                                    iloradius3 = 1; iupradius3 = 1; jloradius3 = 1; jupradius3 = 1;

                                    while (dhobliquemax2 > 0 && available_for_delta_m > 0)
                                    // starting in the lowest cell, we will now find the lowest higher nb
                                    // we will continue looking for lowest higher nbs as long as we find one and have sediment left
                                    // while doing this, we will raise our delta with us
                                    {
                                        find_lowest_higher_oblique_neighbour(active_depression);

                                        if (dhobliquemax2 < 0)
                                        {   // if we have found a lower neighbour
                                            // leave the current delta and bring the remaining sediment to the lower oblique neighbour
                                            if (diagnostic_mode == 1) { Debug.WriteLine(" lowest oblique neighbour is lower - moving sediment down "); }
                                            cleardelta(iloradius3, iupradius3, jloradius3, jupradius3, rowlowestobnb, collowestobnb);
                                            if (iloradius2 < -(rowlowestobnb + II - startrow)) { iloradius2 = -(rowlowestobnb + II - startrow); }
                                            if (iupradius2 < (rowlowestobnb + II - startrow)) { iupradius2 = (rowlowestobnb + II - startrow); }
                                            if (jloradius2 < -(collowestobnb + JJ - startcol)) { jloradius2 = -(collowestobnb + JJ - startcol); }
                                            if (jupradius2 < (collowestobnb + JJ - startcol)) { jupradius2 = (collowestobnb + JJ - startcol); }
                                            dhobliquemax1 = ((GlobalMethods.dtm[startrow, startcol] + GlobalMethods.dz_ero_m[startrow, startcol] + +GlobalMethods.dz_sed_m[startrow, startcol]) - (GlobalMethods.dtm[rowlowestobnb + II, collowestobnb + JJ] + GlobalMethods.dz_ero_m[rowlowestobnb + II, collowestobnb + JJ] + GlobalMethods.dz_sed_m[rowlowestobnb + II, collowestobnb + JJ]) - (Math.Sqrt((startrow - rowlowestobnb - II) * (startrow - rowlowestobnb - II) + (startcol - collowestobnb - JJ) * (startcol - collowestobnb - JJ)) * GlobalMethods.dx * tangent_of_delta)) - 0.0000001;
                                        }
                                        if (dhobliquemax2 == 0)
                                        {
                                            // this is not supposed to happen because in this case the search in find_lowest_higher_oblique_nb must go on.
                                            if (diagnostic_mode == 1) { Debug.WriteLine("Warning. Found dhobliquemax = 0 outside of find_lowest_higher_ob_nb"); }
                                        }
                                        if (dhobliquemax2 > 0)
                                        {
                                            if (diagnostic_mode == 1) { Debug.WriteLine(" lowest oblique neighbour is higher - raising delta"); }
                                            if (diagnostic_mode == 1) { Debug.WriteLine(" available " + available_for_delta_m + "m and space for " + (deltasize * dhobliquemax2) + " m"); }
                                            if (available_for_delta_m >= deltasize * dhobliquemax2)
                                            {
                                                if (diagnostic_mode == 1) { Debug.WriteLine(" raising delta to higher oblique level "); }
                                                raise_delta_completely(active_depression);
                                            }
                                            if (available_for_delta_m < deltasize * dhobliquemax2)
                                            {
                                                if (diagnostic_mode == 1) { Debug.WriteLine(" raising delta as far as possible given sediment "); }
                                                raise_delta_partly(active_depression);
                                                if (diagnostic_mode == 1) { minimaps(row, col); }
                                                if (obnbchanged == 0) { cleardelta(iloradius3, iupradius3, jloradius3, jupradius3, rowlowestobnb, collowestobnb); }
                                                // if the starting cell was raised above lakelevel, it is no longer member of the lake, and we have taken care of that in raise_delta_partly by changing obnb. 
                                                // this must not be removed, so if obnbchanged != 0, we do not clear the delta.
                                            }

                                        } // end if dhoblmax2 > 0
                                    } // end while dhobliquemax2 > 0
                                } // end while sediment_available
                                cleardelta(iloradius3, iupradius3, jloradius3, jupradius3, rowlowestobnb, collowestobnb);
                            } //end if not outlet 
                        } // end if depressio
                    } // end if boundaries
                } // end for col
            } //end for row 

            // we now divide the total amount of extra water (the amount replaced by sediment in the lake) over the (max 5) outlets.
            int outletcounter = 0;
            while (GlobalMethods.drainingoutlet_col[active_depression, outletcounter] != -1)
            {
                outletcounter++;
                if (outletcounter == 5) { break; }
            }
            for (i = 0; i < outletcounter; i++)
            {
                GlobalMethods.waterflow_m3[GlobalMethods.drainingoutlet_col[active_depression, i], GlobalMethods.drainingoutlet_col[active_depression, i]] += GlobalMethods.dx * GlobalMethods.dx * (depressionsum_sediment_m) / outletcounter;
            }
            diagnostic_mode = 1;
        }

        void find_lowest_oblique_neighbour(int this_depression) // to determine where to start or continue with current delta 
        {
            if (GlobalMethods.t > 300000) { diagnostic_mode = 1; }
            if (diagnostic_mode == 1) { Debug.WriteLine(" entered find_lowest_oblique_neighbour"); }
            // finds the lowest oblique neighbour of the current delta
            // affects (changes) global doubles dhobliquemax1 et al
            //int this_depression = Math.Abs(depression[startrow, startcol]);

            int readysearching = 0;
            while (readysearching == 0)
            {
                readysearching = 1;      // we expect to be ready searching, but when not, we will set this to 0
                if (diagnostic_mode == 1) { Debug.WriteLine(" dhobliquemax1 : " + dhobliquemax1); }
                if (diagnostic_mode == 1) { Debug.WriteLine(" ilo " + iloradius2 + ", iup " + iupradius2 + ", jlo " + jloradius2 + ", jup " + jupradius2 + ", row: " + startrow + ", col " + startcol); }
                for (i = -1 * iloradius2; i <= iupradius2; i++)
                {
                    for (j = -1 * jloradius2; j <= jupradius2; j++)
                    {
                        if ((startrow + i >= 0) && (startrow + i < GlobalMethods.nr) && (startcol + j >= 0) && (startcol + j < GlobalMethods.nc) && !((i == 0) && (j == 0)) && GlobalMethods.dtm[startrow + i, startcol + j] != -9999) //&& !((startrow + i == row) && (startcol + j == col))
                        { // boundary check while looking around startrow startcol for the neighbours of the entire current delta
                            if (diagnostic_mode == 1) { Debug.WriteLine(" oblique neighbour now " + (startrow + i) + " " + (startcol + j) + ", depression " + depression[startrow + i, startcol + j] + " GlobalMethods.dtm " + (GlobalMethods.dtm[startrow + i, startcol + j] + GlobalMethods.dz_ero_m[startrow + i, startcol + j] + GlobalMethods.dz_sed_m[startrow + i, startcol + j])); }
                            if (depression[startrow + i, startcol + j] == this_depression || depression[startrow + i, startcol + j] == -this_depression)
                            {  // if you are a member of this depression
                                dhoblique = (GlobalMethods.dtm[startrow, startcol] + GlobalMethods.dz_ero_m[startrow, startcol] + GlobalMethods.dz_sed_m[startrow, startcol]) - (GlobalMethods.dtm[startrow + i, startcol + j] + GlobalMethods.dz_ero_m[startrow + i, startcol + j] + GlobalMethods.dz_sed_m[startrow + i, startcol + j]) - (Math.Sqrt((i * i) + (j * j)) * GlobalMethods.dx * tangent_of_delta);
                                if (diagnostic_mode == 1) { Debug.WriteLine(" oblique neighbour now " + (startrow + i) + " " + (startcol + j) + " :" + dhoblique); }
                                //if (diagnostic_mode == 1) { Debug.WriteLine(" dhobliquemax1 : " + dhobliquemax1); }
                                if ((dhoblique > dhobliquemax1))
                                {      // vanwege afkortinsverschillen 0.000000000 etc 1
                                    dhobliquemax1 = dhoblique;
                                    rowlowestobnb = startrow + i;
                                    collowestobnb = startcol + j;
                                    //deltasize = 1;
                                    if (this_depression > 0)
                                    {
                                        if (diagnostic_mode == 1) { Debug.WriteLine(" lowest oblique neighbour now " + rowlowestobnb + " " + collowestobnb + " dhobliquemax1 : " + dhobliquemax1); }
                                    }
                                    lower_nb_exists = 1;
                                    readysearching = 0;
                                    if (i == -1 * iloradius2) { iloradius2++; }
                                    if (i == iupradius2) { iupradius2++; }
                                    if (j == -1 * jloradius2) { jloradius2++; }
                                    if (j == jupradius2) { jupradius2++; }
                                    if (diagnostic_mode == 1) { Debug.WriteLine(" ilo " + iloradius2 + ", iup " + iupradius2 + ", jlo " + jloradius2 + ", ju2 " + jupradius2); }
                                } // end if dhoblique  < dhobliquemax1)
                                if (dhoblique < 0.0000000001 && dhobliquemax1 < 0.0000000001 && dhoblique > -0.0000000001)
                                { 	// in this case, we may have filled the present neighbour
                                    // in an earlier stage, but have had to clear the delta because
                                    // a lower obnb was found. We must look beyond this equally-high
                                    // oblique nb to be able to find this lower obnb.....
                                    if (diagnostic_mode == 1) { Debug.WriteLine("found previous delta"); }
                                    if (i == -1 * iloradius2) { iloradius2++; readysearching = 0; }
                                    if (i == iupradius2) { iupradius2++; readysearching = 0; }
                                    if (j == -1 * jloradius2) { jloradius2++; readysearching = 0; }
                                    if (j == jupradius2) { jupradius2++; readysearching = 0; }
                                } // end if dhoblique and dhobliquemax1 are zero
                            } // end if depression = depression
                        } // end if boundaries
                    } // end for j
                } // end if
            } // end while readysearching = 0
            if (diagnostic_mode == 1) { Debug.WriteLine(" ready searching - dhobliquemax1: " + dhobliquemax1 + ", GlobalMethods.dx: " + GlobalMethods.dx + ", row " + rowlowestobnb + ", col " + collowestobnb); }
        }

        void find_lowest_higher_oblique_neighbour(int here_depression) // to determine to which level (and cell) the current delta can be raised 
        {
            if (GlobalMethods.t > 1000000) { diagnostic_mode = 1; }
            if (diagnostic_mode == 1) { Debug.WriteLine(" entered find_lowest_higher_oblique_neighbour"); }
            if (diagnostic_mode == 1) { Debug.WriteLine(" rowlow = " + rowlowestobnb + " collow = " + collowestobnb + " range " + iloradius3 + iupradius3 + jloradius3 + jupradius3); }
            if (diagnostic_mode == 1 && rowlowestobnb == 224) { minimaps(rowlowestobnb, collowestobnb); }
            readysearching = 0;
            while (readysearching == 0)
            {
                dhobliquemax2 = 99999.99;
                for (i = -1 * iloradius3; i <= iupradius3; i++)
                {
                    for (j = -1 * jloradius3; j <= jupradius3; j++)
                    {
                        if (((rowlowestobnb + i) >= 0) && ((rowlowestobnb + i) < GlobalMethods.nr) && ((collowestobnb + j) >= 0) && ((collowestobnb + j) < GlobalMethods.nc) && !((i == 0) && (j == 0)))
                        { // boundaries
                            if (depression[rowlowestobnb + i, collowestobnb + j] == here_depression)
                            {
                                dhoblique = -((GlobalMethods.dtm[rowlowestobnb, collowestobnb]) + GlobalMethods.dz_ero_m[rowlowestobnb, collowestobnb] + GlobalMethods.dz_sed_m[rowlowestobnb, collowestobnb]) + (GlobalMethods.dtm[rowlowestobnb + i, collowestobnb + j] + GlobalMethods.dz_ero_m[rowlowestobnb + i, collowestobnb + j] + GlobalMethods.dz_sed_m[rowlowestobnb + i, collowestobnb + j]) + (Math.Sqrt(Math.Pow((rowlowestobnb + i - startrow), 2) + Math.Pow((collowestobnb + j - startcol), 2)) - Math.Sqrt(Math.Pow((rowlowestobnb - startrow), 2) + Math.Pow((collowestobnb - startcol), 2))) * GlobalMethods.dx * tangent_of_delta;
                                if (dhoblique != 0 && dhoblique < dhobliquemax2)
                                {
                                    readysearching = 1;
                                    dhobliquemax2 = dhoblique;
                                    II = i; JJ = j;
                                    if (diagnostic_mode == 1 && dhobliquemax2 < 0) { Debug.WriteLine("cell " + (rowlowestobnb + i) + " " + (collowestobnb + j) + " is a lower obnb (dhoblmax2 = " + dhobliquemax2 + ") from " + rowlowestobnb + collowestobnb); }
                                    // 
                                } // end if
                                if (dhoblique < 0.0000000001 && dhoblique > -0.0000000001)
                                { 		//this means we have transferred to a 'new delta' after encountering a negative dhobliquemax2 first,
                                    //then filled this new delta to the point that we are as high as the previous one
                                    // we therefore encounter one or more cells with dhoblique = 0
                                    //we must incorporate these cells into the delta and increase the searchradius
                                    depression[rowlowestobnb + i, collowestobnb + j] = -here_depression;
                                    deltasize++;
                                    readysearching = 0;
                                    if (diagnostic_mode == 1) { Debug.WriteLine("added " + (rowlowestobnb + i) + " " + (collowestobnb + j) + " (dhoblique 0.0000) to the delta and will increase searchradius for lowest higher nbour"); }
                                } //end if dhoblmax2 == 0
                            } // end if depression == depression
                        } // end if bndries
                    } // end for j
                } // end for i
                if (readysearching == 0)
                {
                    // If no neighbour within the original radius3 from obnb was member of a lake (which can happen if they were raised to dtm_fill earlier), then 
                    // we will increase the searchradius (there can be another lake cell somewhere, because the lake is not being filled completely).
                    iloradius3++; if (diagnostic_mode == 1) { Debug.WriteLine("ilo is higher"); }
                    iupradius3++; if (diagnostic_mode == 1) { Debug.WriteLine("iup is higher"); }
                    jloradius3++; if (diagnostic_mode == 1) { Debug.WriteLine("jlo is higher"); }
                    jupradius3++; if (diagnostic_mode == 1) { Debug.WriteLine("jup is higher"); }
                    //however, no higher ob nb is necessarily present. The current obnb, and possibly its delta, need not have a lake-neighbour with a higher dhoblique
                    //this is for instance the case when we are currently looking at the last cell in the lake.
                    //In that case, we must fill the current delta as much as possible with the existing sediment. A corresponding dhoblique must be sent back to control.
                    if (rowlowestobnb - iloradius3 <= iloedge[here_depression] && rowlowestobnb + iupradius3 >= iupedge[here_depression] &&
                        collowestobnb - jloradius3 <= jloedge[here_depression] && collowestobnb + jupradius3 >= jupedge[here_depression])
                    {
                        readysearching = 1;
                        dhobliquemax2 = (GlobalMethods.dtmfill_A[rowlowestobnb, collowestobnb] - (GlobalMethods.dtm[rowlowestobnb, collowestobnb] + GlobalMethods.dz_ero_m[rowlowestobnb, collowestobnb] + GlobalMethods.dz_sed_m[rowlowestobnb, collowestobnb]));
                        if (dhobliquemax2 == 0) { dhobliquemax2 = 1; }//ArT If the dhoblique is zero, because the obnb was the outlet with dtmfill== GlobalMethods.dtm, then we give an emergency value to dhobliquemax2  
                        if (diagnostic_mode == 1) { Debug.WriteLine(" search for lower ob nb finished - not found - sending dhobliquemax as equal to fillspace of " + dhobliquemax2); }
                    }
                }
            } //end while readysearching == 1
        }

        void raise_delta_completely(int this_depression) // raise delta completely (and then go on raising it to higher obl heights)  
        {
            int size;
            // the amount of sed_in_trans left is enough to raise the entire delta with dhobliquemax1
            if (diagnostic_mode == 1) { Debug.WriteLine(" 2: to be added: " + dhobliquemax2 + ", sed_for_delta: " + available_for_delta_m + "m, deltasize: " + deltasize); }
            if (diagnostic_mode == 1) { Debug.WriteLine(" rowlowestobnb " + rowlowestobnb + " collowestobnb " + collowestobnb + " and " + iloradius3 + iupradius3 + jloradius3 + jupradius3); }
            for (i = -1 * iloradius3; i <= iupradius3; i++)
            {
                for (j = -1 * jloradius3; j <= jupradius3; j++)
                {
                    if (((rowlowestobnb + i) >= 0) && ((rowlowestobnb + i) < GlobalMethods.nr) && ((collowestobnb + j) >= 0) && ((collowestobnb + j) < GlobalMethods.nc) && !(rowlowestobnb + i == row && collowestobnb + j == col) && GlobalMethods.dtm[rowlowestobnb + i, collowestobnb + j] != -9999)
                    { // boundaries, note that i==0 && j==0 is allowed  ;we can raise rowlowestobnb,colloewsobnb when it is part of the delta.
                        if (depression[rowlowestobnb + i, collowestobnb + j] == -this_depression) // i.e. if cell is part of present delta
                        {
                            GlobalMethods.dtm[rowlowestobnb + i, collowestobnb + j] += dhobliquemax2;
                            available_for_delta_m -= dhobliquemax2;
                            if (GlobalMethods.dtm[rowlowestobnb + i, collowestobnb + j] == -1) { Debug.WriteLine("A1 cell " + (rowlowestobnb + i) + " " + (collowestobnb + j) + " has an altitude of -1 now"); minimaps((rowlowestobnb + i), (collowestobnb + j)); }
                            if (available_for_delta_m < 0) { Debug.WriteLine(" Error: negative sediment for delta " + available_for_delta_m + " m"); minimaps(rowlowestobnb + i, collowestobnb + j); }
                            if (diagnostic_mode == 1) { Debug.WriteLine(" raised " + (rowlowestobnb + i) + " " + (collowestobnb + j) + " with " + dhobliquemax2 + " to " + (GlobalMethods.dtm[rowlowestobnb + i, collowestobnb + j] + GlobalMethods.dz_ero_m[rowlowestobnb + i, collowestobnb + j] + GlobalMethods.dz_sed_m[rowlowestobnb + i, collowestobnb + j]) + " sed_for_delta now " + available_for_delta_m); }
                            if (diagnostic_mode == 1 && dhobliquemax2 < 0) { Debug.WriteLine(" warning: dhobliquemax2 is less than zero at " + (rowlowestobnb + i) + " " + (collowestobnb + j) + " " + dhobliquemax2); }
                            //if (diagnostic_mode == 1) { MessageBox.Show("warning - extremely high coarse sed_in_trans:" + GlobalMethods.sediment_in_transport_kg[startrow, startcol,0]); }
                            if ((GlobalMethods.dtm[rowlowestobnb + i, collowestobnb + j] + GlobalMethods.dz_ero_m[rowlowestobnb + i, collowestobnb + j] + GlobalMethods.dz_sed_m[rowlowestobnb + i, collowestobnb + j]) > GlobalMethods.dtmfill_A[rowlowestobnb + i, collowestobnb + j])
                            {   // then we have raised this cell too high
                                if (diagnostic_mode == 1) { Debug.WriteLine("1 we change the altitude of " + (rowlowestobnb + i) + " " + (collowestobnb + j) + " (depressionlevel " + depressionlevel[this_depression] + ") from " + GlobalMethods.dtm[rowlowestobnb + i, collowestobnb + j] + " to " + GlobalMethods.dtmfill_A[rowlowestobnb + i, collowestobnb + j]); }
                                available_for_delta_m += ((GlobalMethods.dtm[rowlowestobnb + i, collowestobnb + j] + GlobalMethods.dz_ero_m[rowlowestobnb + i, collowestobnb + j] + GlobalMethods.dz_sed_m[rowlowestobnb + i, collowestobnb + j]) - GlobalMethods.dtmfill_A[rowlowestobnb + i, collowestobnb + j]);
                                double[] local_s_i_t_kg = new double[5] { 0, 0, 0, 0, 0 };
                                for (size = 0; size < GlobalMethods.n_texture_classes; size++)
                                {
                                    //any sediment in transport that was possibly waiting for consideration later than the current startrow, startcol is taken into this startrow startcol to make a bigger delta here
                                    local_s_i_t_kg[size] = GlobalMethods.sediment_in_transport_kg[rowlowestobnb + i, collowestobnb + j, size];
                                    GlobalMethods.sediment_in_transport_kg[rowlowestobnb + i, collowestobnb + j, size] = 0;
                                }
                                if (diagnostic_mode == 1) { Debug.WriteLine(" adding " + calc_thickness_from_mass(local_s_i_t_kg, 0, 0) + " to delta-available " + available_for_delta_m); }
                                available_for_delta_m += calc_thickness_from_mass(local_s_i_t_kg, 0, 0);
                                if (available_for_delta_m < 0) { Debug.WriteLine("5 negative sediment for delta " + available_for_delta_m + " m"); }
                                GlobalMethods.lake_sed_m[rowlowestobnb + i, collowestobnb + j] += ((GlobalMethods.dtm[rowlowestobnb + i, collowestobnb + j] + GlobalMethods.dz_ero_m[rowlowestobnb + i, collowestobnb + j] + GlobalMethods.dz_sed_m[rowlowestobnb + i, collowestobnb + j]) - GlobalMethods.dtmfill_A[rowlowestobnb + i, collowestobnb + j]);
                                if (GlobalMethods.lake_sed_m[rowlowestobnb + i, collowestobnb + j] < -0.0000001) { Debug.WriteLine("1 Warning: negative lake deposition in " + (rowlowestobnb + i) + " " + (collowestobnb + j) + " of " + GlobalMethods.lake_sed_m[rowlowestobnb + i, collowestobnb + j] + " GlobalMethods.dtm " + GlobalMethods.dtm[rowlowestobnb + i, collowestobnb + j] + " fill " + GlobalMethods.dtmfill_A[rowlowestobnb + i, collowestobnb + j]); minimaps(rowlowestobnb + i, collowestobnb + j); }
                                GlobalMethods.dtm[rowlowestobnb + i, collowestobnb + j] = GlobalMethods.dtmfill_A[rowlowestobnb + i, collowestobnb + j] - GlobalMethods.dz_ero_m[rowlowestobnb + i, collowestobnb + j] - GlobalMethods.dz_sed_m[rowlowestobnb + i, collowestobnb + j];
                                if (GlobalMethods.dtm[rowlowestobnb + i, collowestobnb + j] == -1) { Debug.WriteLine("A2 cell " + (rowlowestobnb + i) + " " + (collowestobnb + j) + " has an altitude of -1 now"); minimaps((rowlowestobnb + i), (collowestobnb + j)); }
                                depression[rowlowestobnb + i, collowestobnb + j] = 0;
                                deltasize--;
                                if (diagnostic_mode == 1) { Debug.WriteLine("1: " + (rowlowestobnb + i) + " " + (collowestobnb + j) + " (depressionlevel " + depressionlevel[depression[row, col]] + ") now at " + GlobalMethods.dtm[rowlowestobnb + i, collowestobnb + j] + " = fill_A " + GlobalMethods.dtmfill_A[rowlowestobnb + i, collowestobnb + j] + " sed for delta " + available_for_delta_m); }
                                if (diagnostic_mode == 1) { Debug.WriteLine("decreased deltasize with 1 to " + deltasize); }
                            } // end if GlobalMethods.dtm > depressionlevel
                        }   // end if member of delta
                    }  // end if boundaries
                }  // end for j
            }  // end for i
            if (diagnostic_mode == 1)
            {
                Debug.WriteLine(" raised and now move on to enlarge delta with " + (rowlowestobnb + II) + " " + (collowestobnb + JJ) + " available for delta now " + available_for_delta_m + "m");
                Debug.WriteLine(" by changing " + (rowlowestobnb + II) + "," + (collowestobnb + JJ) + " from " + depression[rowlowestobnb + II, collowestobnb + JJ] + " into " + (-this_depression));
            }
            depression[rowlowestobnb + II, collowestobnb + JJ] = -this_depression;
            deltasize++;
            if (diagnostic_mode == 1) { Debug.WriteLine(" raised deltasize with 1 to " + deltasize + ". " + (rowlowestobnb + II) + " " + (collowestobnb + JJ) + " now negative lakevalue"); }
            //dhobliquemax2 = 0;
            if (II == -iloradius3) { iloradius3++; }
            if (II == iupradius3) { iupradius3++; }
            if (JJ == -jloradius3) { jloradius3++; }
            if (JJ == jupradius3) { jupradius3++; }
        }

        void raise_delta_partly(int this_depression) // raise partly (and then go to new cells on border of lake)  
        {

            int size;
            //Debug.WriteLine( "raising delta partly for dep " + this_depression);
            mem_m = available_for_delta_m / deltasize;
            available_for_delta_m = 0;
            if (diagnostic_mode == 1) { Debug.WriteLine(" 1: to be added: " + mem_m + ", deltasize: " + deltasize); }
            tempx = rowlowestobnb; tempy = collowestobnb;
            obnbchanged = 0;
            for (i = -1 * iloradius3; i <= iupradius3; i++)
            {
                for (j = -1 * jloradius3; j <= jupradius3; j++)
                {
                    if (((tempx + i) >= 0) && ((tempx + i) < GlobalMethods.nr) && ((tempy + j) >= 0) && ((tempy + j) < GlobalMethods.nc) && !(tempx + i == row && tempy + j == col) && GlobalMethods.dtm[tempx + i, tempy + j] != -9999)
                    { // boundaries
                        if (depression[tempx + i, tempy + j] == -this_depression)
                        {
                            GlobalMethods.dtm[tempx + i, tempy + j] += mem_m;
                            GlobalMethods.lake_sed_m[tempx + i, tempy + j] += mem_m;
                            if (GlobalMethods.dtm[rowlowestobnb + i, collowestobnb + j] == -1) { Debug.WriteLine("B cell " + (rowlowestobnb + i) + " " + (collowestobnb + j) + " has an altitude of -1 now"); minimaps((rowlowestobnb + i), (collowestobnb + j)); }
                            if (GlobalMethods.lake_sed_m[tempx + i, tempy + j] < -0.0000001) { Debug.WriteLine("4 Warning: negative lake deposition in " + (tempx + i) + " " + (tempy + j) + " of " + GlobalMethods.lake_sed_m[tempx + i, tempy + j]); minimaps(tempx + i, tempy + j); }
                            if (diagnostic_mode == 1) { Debug.WriteLine(" added " + mem_m + " to cell " + (tempx + i) + " " + (tempy + j)); }
                            if ((GlobalMethods.dtm[tempx + i, tempy + j] + GlobalMethods.dz_ero_m[tempx + i, tempy + j] + GlobalMethods.dz_sed_m[tempx + i, tempy + j]) > GlobalMethods.dtmfill_A[tempx + i, tempy + j])
                            {
                                if (diagnostic_mode == 1) { Debug.WriteLine(" cell " + (tempx + i) + " " + (tempy + j) + " raised above filllevel " + GlobalMethods.dtmfill_A[tempx + i, tempy + j] + ", to " + (GlobalMethods.dtm[tempx + i, tempy + j] + GlobalMethods.dz_ero_m[tempx + i, tempy + j] + GlobalMethods.dz_sed_m[tempx + i, tempy + j])); }
                                available_for_delta_m += ((GlobalMethods.dtm[tempx + i, tempy + j] + GlobalMethods.dz_ero_m[tempx + i, tempy + j] + GlobalMethods.dz_sed_m[tempx + i, tempy + j]) - GlobalMethods.dtmfill_A[tempx + i, tempy + j]);
                                for (size = 0; size < GlobalMethods.n_texture_classes; size++)
                                {
                                    local_s_i_t_kg[size] = GlobalMethods.sediment_in_transport_kg[tempx + i, tempy + j, size];
                                    GlobalMethods.sediment_in_transport_kg[tempx + i, tempy + j, size] = 0;   //ART recently changed, should solve a bug (this line did not exist, violation mass balance
                                }
                                available_for_delta_m += calc_thickness_from_mass(local_s_i_t_kg, 0, 0);
                                if (available_for_delta_m < 0) { Debug.WriteLine("9 negative sediment in transport (m) remaining for delta " + available_for_delta_m + "m"); }
                                if (diagnostic_mode == 1) { Debug.WriteLine(" A we change the altitude of " + (tempx + i) + " " + (tempy + j) + " (depressionlevel " + depressionlevel[this_depression] + ") from " + (GlobalMethods.dtm[tempx + i, tempy + j] + GlobalMethods.dz_ero_m[tempx + i, tempy + j] + GlobalMethods.dz_sed_m[tempx + i, tempy + j]) + " to " + GlobalMethods.dtmfill_A[tempx + i, tempy + j]); }
                                if (tempx + i == row && tempy + j == col) { Debug.WriteLine("we are changing outlet " + tempx + " " + tempy + " into 0"); }
                                GlobalMethods.lake_sed_m[tempx + i, tempy + j] -= ((GlobalMethods.dtm[tempx + i, tempy + j] + GlobalMethods.dz_ero_m[tempx + i, tempy + j] + GlobalMethods.dz_sed_m[tempx + i, tempy + j]) - GlobalMethods.dtmfill_A[tempx + i, tempy + j]);
                                if (GlobalMethods.lake_sed_m[tempx + i, tempy + j] < -0.0000001) { Debug.WriteLine("3 Warning: negative lake deposition in " + (tempx + i) + " " + (tempy + j) + " of " + GlobalMethods.lake_sed_m[tempx + i, tempy + j] + " alt " + (GlobalMethods.dtm[tempx + i, tempy + j] + GlobalMethods.dz_ero_m[tempx + i, tempy + j] + GlobalMethods.dz_sed_m[tempx + i, tempy + j]) + " fill " + GlobalMethods.dtmfill_A[tempx + i, tempy + j]); minimaps(tempx + i, tempy + j); }
                                GlobalMethods.dtm[tempx + i, tempy + j] = (GlobalMethods.dtmfill_A[tempx + i, tempy + j] - GlobalMethods.dz_ero_m[tempx + i, tempy + j] - GlobalMethods.dz_sed_m[tempx + i, tempy + j]); //so that with ero and sed, it equals dtmfill
                                if (GlobalMethods.dtm[tempx + i, tempy + j] == -1) { Debug.WriteLine("C cell " + (tempx + i) + " " + (tempy + j) + " has an altitude of -1 now"); minimaps(tempx + i, tempy + j); } //
                                if (diagnostic_mode == 1) { Debug.WriteLine(" will change depressionmembership of " + (tempx + i) + " " + (tempy + j) + " from " + depression[tempx + i, tempy + j] + " to 0"); }
                                if (diagnostic_mode == 1) { Debug.WriteLine(" II = " + II + ", JJ = " + JJ); }
                                depression[tempx + i, tempy + j] = 0;
                                if (diagnostic_mode == 1) { Debug.WriteLine(" obnbchanged? - row " + row + " col " + col + " startrow " + startrow + " startcol " + startcol + " rowobnb " + rowlowestobnb + " colobnb " + collowestobnb + " tempx+i " + (tempx + i) + " tempy+j " + (tempy + j)); }
                                obnbchanged = 1;  // if there is at least one cell that has been raised above dtmfill, then that is the lowest oblique neighbour: the cell with rowlowestobnb, collowestobnb. 
                                // In that case, we must move to a new rowlowestobnb collowestobnb to build the remainder of the delta from there. There is a remainder because some of the sediment used to raise the original
                                // lowest oblique neighbour above dtmfill has been added to available_for_delta again.
                                deltasize--;
                                if (diagnostic_mode == 1) { Debug.WriteLine("decreased deltasize with 1 to " + deltasize); }
                                break;
                            } // end if higher than depressionlevel
                        }  //end if member of delta
                    } // end for boundaries
                } //end for j
            }  // end for i
            if (obnbchanged == 1)
            {
                // this code makes sure that the starting situation for the delta to be built with the remainder of the sediment is correct.
                // this desired situation is: 
                // a) correct amount of sediment (already guaranteed)
                // b) correct cells member of delta (already correct if deltasize > 0, but if deltasize == 0 after removal of the too-high cell, then not correct).
                // c) correct starting cell for delta building (rowlowestobnb and collowestobnb). Not correct because it may be that cell which was just removed from lake and delta. 
                // In case deltasize > 0 and rowlobnb,collobnb is no longer part of the delta, we select any delta-neighbour of the old starting cell. In case deltasize = 0, we select the lowest higher oblique nb for this.

                if (deltasize == 0)
                {
                    deltasize++;
                    if (diagnostic_mode == 1) { Debug.WriteLine("increased deltasize with 1 to " + deltasize); }
                    if (diagnostic_mode == 1) { Debug.WriteLine("lowest oblique neighbour now " + rowlowestobnb + " " + collowestobnb + ", will be: " + (tempx + II) + " " + (tempy + JJ)); }
                    rowlowestobnb = tempx + II;
                    collowestobnb = tempy + JJ;
                    if (diagnostic_mode == 1) { Debug.WriteLine(" will change depressionmembership of " + rowlowestobnb + " " + collowestobnb + " from " + depression[rowlowestobnb, collowestobnb] + " to " + (-this_depression)); }
                    depression[rowlowestobnb, collowestobnb] = -this_depression;
                }
                if (deltasize > 0 && depression[rowlowestobnb, collowestobnb] != -this_depression)
                {
                    if (depression[rowlowestobnb + 1, collowestobnb + 1] == -this_depression) { rowlowestobnb = rowlowestobnb + 1; collowestobnb = collowestobnb + 1; }
                    if (depression[rowlowestobnb + 1, collowestobnb] == -this_depression) { rowlowestobnb = rowlowestobnb + 1; }
                    if (depression[rowlowestobnb + 1, collowestobnb - 1] == -this_depression) { rowlowestobnb = rowlowestobnb + 1; collowestobnb = collowestobnb - 1; }
                    if (depression[rowlowestobnb, collowestobnb + 1] == -this_depression) { collowestobnb = collowestobnb + 1; }
                    if (depression[rowlowestobnb, collowestobnb - 1] == -this_depression) { collowestobnb = collowestobnb - 1; }
                    if (depression[rowlowestobnb - 1, collowestobnb + 1] == -this_depression) { rowlowestobnb = rowlowestobnb - 1; collowestobnb = collowestobnb + 1; }
                    if (depression[rowlowestobnb - 1, collowestobnb] == -this_depression) { rowlowestobnb = rowlowestobnb - 1; }
                    if (depression[rowlowestobnb - 1, collowestobnb - 1] == -this_depression) { rowlowestobnb = rowlowestobnb - 1; collowestobnb = collowestobnb - 1; }
                }
            }
            if (diagnostic_mode == 1) { Debug.WriteLine(" sed_for_delta is now " + available_for_delta_m + " and deltasize = " + deltasize); }
            //diagnostic_mode = 0;
        }

        #endregion

        #region initialisation code

        void initialise_once()        //fills the inputgrids with values
        {

            if (guiVariables.Daily_water && GlobalMethods.input_data_error == false)
            {
                try
                {
                    filename = this.dailyP.Text;
                    if (memory_records_d == false) { makedailyrecords(filename); }
                    // Debug.WriteLine("guiVariables.P_all record created successfully");

                    int namelength = filename.Length;
                    string fn = filename.Substring(0, (namelength - 4));
                    filename = fn + "_sc" + P_scen + ".csv";

                    GlobalMethods.read_record(filename, guiVariables.P_all);
                    // Debug.WriteLine("guiVariables.P_all read successfully");

                    //filename = this.dailyET0.Text;
                    //if (memory_records_d == false) { makedailyrecords(filename); }
                    //read_record(filename, guiVariables.ET0_all);
                    //Debug.WriteLine("guiVariables.ET0_all read successfully");

                    filename = this.dailyD.Text;
                    if (memory_records_d == false) { makedailyrecords(filename); }
                    GlobalMethods.read_record(filename, guiVariables.D_all);
                    // Debug.WriteLine("guiVariables.D_all read successfully");

                    filename = this.dailyT_avg.Text;
                    if (memory_records_d == false) { makedailyrecords(filename); }
                    GlobalMethods.read_record(filename, guiVariables.Tavg_all);
                    // Debug.WriteLine("guiVariables.Tavg_all read successfully");

                    filename = this.dailyT_min.Text;
                    if (memory_records_d == false) { makedailyrecords(filename); }
                    GlobalMethods.read_record(filename, guiVariables.Tmin_all);
                    // Debug.WriteLine("guiVariables.Tmin_all read successfully");

                    filename = this.dailyT_max.Text;
                    if (memory_records_d == false) { makedailyrecords(filename); }
                    GlobalMethods.read_record(filename, guiVariables.Tmax_all);
                    // Debug.WriteLine("guiVariables.Tmax_all read successfully");

                    //// Hargreaves extraterrestrial radiation
                    //// http://www.fao.org/docrep/X0490E/x0490e07.htm#radiation
                    //double dr, delta, ws;
                    //double lat_st = Math.PI/180*(System.Convert.ToDouble(latitude_deg.Text)+ System.Convert.ToDouble(latitude_min.Text) / 60); // latitude in radians

                    //for (double day_ra = 1; day_ra <= 365; day_ra++)
                    //{
                    //    dr = 1 + 0.033 * Math.Cos(2 * Math.PI * (day_ra / 365)); //inverse relative distance Earth-Sun
                    //    delta = 0.409 * Math.Sin(2 * Math.PI * (day_ra / 365) - 1.39); //solar declination [rad].
                    //    ws = Math.Acos(-Math.Tan(lat_st) * Math.Tan(delta)); // sunset hour angle [rad]
                    //    Ra_ann[Convert.ToInt16(day_ra - 1)] = (24 * 60 / Math.PI * 0.082 * dr * (ws * Math.Sin(lat_st) * Math.Sin(delta) + Math.Cos(lat_st) * Math.Cos(delta) * Math.Sin(ws))) * 0.408; // extraterrestrial radiation [mm d-1]
                    //}

                    Ra_rcm = new double[GlobalMethods.nr, GlobalMethods.nc, 12];
                    guiVariables.OFy_m = new double[GlobalMethods.nr, GlobalMethods.nc, 10]; // 0: outflow, 1:8 flow to neighbours, 9: inflow
                    Iy = new double[GlobalMethods.nr, GlobalMethods.nc];
                    waterfactor = new double[GlobalMethods.nr, GlobalMethods.nc];
                    pond_y = new double[GlobalMethods.nr, GlobalMethods.nc];
                    outflow_y = new double[GlobalMethods.nr, GlobalMethods.nc];
                    total_outflow_y = new double[GlobalMethods.nr, GlobalMethods.nc];
                    water_balance_m = new double[GlobalMethods.nr, GlobalMethods.nc, 5]; // 1: rainfall, 2: actual ET, 3: runon, 4: runoff, 5: I
                    ETay = new double[GlobalMethods.nr, GlobalMethods.nc];
                    ET0y = new double[GlobalMethods.nr, GlobalMethods.nc];
                    veg_correction_factor = new double[GlobalMethods.nr, GlobalMethods.nc];
                    GlobalMethods.waterflow_m3 = new double[GlobalMethods.nr, GlobalMethods.nc];
                    vegetation_type = new int[GlobalMethods.nr, GlobalMethods.nc];

                    for (int row = 0; row < GlobalMethods.nr; row++)
                    {
                        for (int col = 0; col < GlobalMethods.nc; col++)
                        {
                            veg_correction_factor[row, col] = 1;
                            vegetation_type[row, col] = 0;
                        }
                    }

                    Ks_topsoil_mh = new double[GlobalMethods.nr, GlobalMethods.nc];
                    Ks_md = new double[GlobalMethods.nr, GlobalMethods.nc, GlobalMethods.max_soil_layers];
                    stagdepth = new double[GlobalMethods.nr, GlobalMethods.nc];
                    snow_m = 0;
                    snow_start_m = 0;
                    snowfall_m = 0;
                    snowmelt_factor_mTd = Convert.ToDouble(snowmelt_factor_textbox.Text);
                    snow_threshold_C = Convert.ToDouble(snow_threshold_textbox.Text);
                    // snowmelt factor now assumed to be 0.004 m per degree per day, as a mean from the following paper: Hock 2003, https://www.sciencedirect.com/science/article/pii/S0022169403002579#BIB50
                    // DEVELOP I didn'GlobalMethods.t take the correct parameter (DDF instead of Fm, see paper). Change to right parameter based on the paper of Gottlieb, 1980. Other developments are varying melt factors, based on month, incoming radiation etc. (See Hock 2003)
                }
                catch { Debug.WriteLine(" problem preparing for hourly water balance "); }
            }

            if (check_space_soildepth.Checked && GlobalMethods.input_data_error == false)
            {
                filename = this.soildepth_input_filename_textbox.Text;
                GlobalMethods.read_double(filename, GlobalMethods.soildepth_m);
                Debug.WriteLine("read soildepth");
            }
            if (check_space_till_fields.Checked && GlobalMethods.input_data_error == false)
            {
                filename = guivariables.Tillfields_input_filename_textbox;
                GlobalMethods.read_integer(filename, GlobalMethods.tillfields);
                Debug.WriteLine("read GlobalMethods.tillfields");
            }
            if (check_space_landuse.Checked && GlobalMethods.input_data_error == false)
            {
                filename = this.landuse_input_filename_textbox.Text;
                GlobalMethods.read_integer(filename, GlobalMethods.landuse);
                Debug.WriteLine("read GlobalMethods.landuse");
            }
            if (check_space_evap.Checked && GlobalMethods.input_data_error == false)
            {
                filename = this.evap_input_filename_textbox.Text;
                GlobalMethods.read_double(filename, GlobalMethods.evapotranspiration);
            }
            if (check_space_infil.Checked && GlobalMethods.input_data_error == false)
            {
                filename = this.infil_input_filename_textbox.Text;
                GlobalMethods.read_double(filename, GlobalMethods.infil);
            }
            if (check_space_rain.Checked && GlobalMethods.input_data_error == false)
            {
                filename = this.rain_input_filename_textbox.Text;
                GlobalMethods.read_double(filename, GlobalMethods.rain);
            }
            // If required, read timeseries instead.
            if (check_time_landuse.Checked && GlobalMethods.input_data_error == false)
            {
                filename = this.landuse_input_filename_textbox.Text;
                GlobalMethods.read_integer(filename, GlobalMethods.landuse);
            }
            if (check_time_evap.Checked && GlobalMethods.input_data_error == false)
            {
                filename = this.evap_input_filename_textbox.Text;
                if (memory_records == false) { makerecords(filename); }
                GlobalMethods.read_record(filename, evap_record);
            }
            if (check_time_infil.Checked && GlobalMethods.input_data_error == false)
            {
                filename = this.infil_input_filename_textbox.Text;
                if (memory_records == false) { makerecords(filename); }
                GlobalMethods.read_record(filename, infil_record);
            }
            if (check_time_rain.Checked && GlobalMethods.input_data_error == false)
            {
                filename = this.rain_input_filename_textbox.Text;
                if (memory_records == false) { makerecords(filename); }
                GlobalMethods.read_record(filename, rainfall_record);
            }
            if (check_time_T.Checked && GlobalMethods.input_data_error == false)
            {
                filename = this.temp_input_filename_textbox.Text;
                if (memory_records == false) { makerecords(filename); }
                GlobalMethods.read_record(filename, temp_record);
            }


            if (guivariables.Check_time_till_fields && GlobalMethods.input_data_error == false) 
            {
                filename = guivariables.Tillfields_input_filename_textbox;
                if (memory_records == false) { makerecords(filename); }
                GlobalMethods.read_record(filename, till_record);
                Debug.WriteLine("Tillage time parameters read");
            }

            try
            {
                // Debug.WriteLine(" assigning starting values for geomorph  ");
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        GlobalMethods.dz_soil[row, col] = 0;
                        if (guiVariables.Creep_Checkbox) { GlobalMethods.sum_creep_grid[row, col] = 0; GlobalMethods.creep[row, col] = 0; }
                        if (guiVariables.Solifluction_checkbox) { GlobalMethods.sum_solifluction[row, col] = 0; }
                        if (guiVariables.Water_ero_checkbox && guiVariables.Only_waterflow_checkbox == false) { GlobalMethods.sum_water_erosion[row, col] = 0; total_sed_export = 0; }
                        if (guiVariables.Treefall_checkbox) { GlobalMethods.dz_treefall[row, col] = 0; GlobalMethods.treefall_count[row, col] = 0; }
                        if (guiVariables.Biological_weathering_checkbox) { GlobalMethods.sum_biological_weathering[row, col] = 0; }
                        if (guiVariables.Frost_weathering_checkbox) { GlobalMethods.sum_frost_weathering[row, col] = 0; }
                        if (guiVariables.Tillage_checkbox) { GlobalMethods.sum_tillage[row, col] = 0; total_sum_tillage = 0; GlobalMethods.dz_till_bd[row, col] = 0; }
                        if (guiVariables.Landslide_checkbox) { GlobalMethods.sum_landsliding[row, col] = 0; total_sum_tillage = 0; }
                        if (GlobalMethods.soildepth_m[row, col] < 0.0 && GlobalMethods.soildepth_m[row, col] != -9999) { soildepth_error += GlobalMethods.soildepth_m[row, col]; GlobalMethods.soildepth_m[row, col] = 0; }
                        if (guiVariables.Ulift_active_checkbox) { GlobalMethods.sum_uplift[row, col] = 0; GlobalMethods.sum_uplift = 0; }
                        if (guiVariables.Tilting_active_checkbox) { GlobalMethods.sum_tilting[row, col] = 0; total_sum_tilting = 0; }
                        if (guiVariables.Check_space_soildepth != true) { GlobalMethods.soildepth_m[row, col] = soildepth_value; }
                        if (guiVariables.Check_space_till_fields != true && guiVariables.Tillage_checkbox)
                        {
                            GlobalMethods.tillfields[row, col] = 1;

                        }



                        if (guiVariables.Water_ero_checkbox && guiVariables.Only_waterflow_checkbox == false)
                        {
                            GlobalMethods.K_fac[row, col] = advection_erodibility; GlobalMethods.P_fac[row, col] = P_act;
                        } //WVG GlobalMethods.K_fac matrix initialisation is needed when GlobalMethods.landuse is disabled

                        if (guiVariables.Water_ero_checkbox && guiVariables.Check_space_landuse == true)
                        {
                            //currently, this will throw an exception if GlobalMethods.landuse is actually spatial //development required //ArT
                            if (GlobalMethods.landuse[row, col] == 1)
                            {
                                GlobalMethods.infil[row, col] *= System.Convert.ToDouble(landuse_determinator.LU1_Inf_textbox.Text);
                                GlobalMethods.evapotranspiration[row, col] *= System.Convert.ToDouble(landuse_determinator.LU1_Evap_textbox.Text);
                                GlobalMethods.K_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU1_Ero_textbox.Text);
                                GlobalMethods.P_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU1_Dep_textbox.Text);
                            }
                            if (GlobalMethods.landuse[row, col] == 2)
                            {
                                GlobalMethods.infil[row, col] *= System.Convert.ToDouble(landuse_determinator.LU2_Inf_textbox.Text);
                                GlobalMethods.evapotranspiration[row, col] *= System.Convert.ToDouble(landuse_determinator.LU2_Evap_textbox.Text);
                                GlobalMethods.K_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU2_Ero_textbox.Text);
                                GlobalMethods.P_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU2_Dep_textbox.Text);
                            }
                            if (GlobalMethods.landuse[row, col] == 3)
                            {
                                GlobalMethods.infil[row, col] *= System.Convert.ToDouble(landuse_determinator.LU3_Inf_textbox.Text);
                                GlobalMethods.evapotranspiration[row, col] *= System.Convert.ToDouble(landuse_determinator.LU3_Evap_textbox.Text);
                                GlobalMethods.K_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU3_Ero_textbox.Text);
                                GlobalMethods.P_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU3_Dep_textbox.Text);
                            }
                            if (GlobalMethods.landuse[row, col] == 4)
                            {
                                GlobalMethods.infil[row, col] *= System.Convert.ToDouble(landuse_determinator.LU4_Inf_textbox.Text);
                                GlobalMethods.evapotranspiration[row, col] *= System.Convert.ToDouble(landuse_determinator.LU4_Evap_textbox.Text);
                                GlobalMethods.K_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU4_Ero_textbox.Text);
                                GlobalMethods.P_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU4_Dep_textbox.Text);
                            }
                            if (GlobalMethods.landuse[row, col] == 5)
                            {
                                GlobalMethods.infil[row, col] *= System.Convert.ToDouble(landuse_determinator.LU5_Inf_textbox.Text);
                                GlobalMethods.evapotranspiration[row, col] *= System.Convert.ToDouble(landuse_determinator.LU5_Evap_textbox.Text);
                                GlobalMethods.K_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU5_Ero_textbox.Text);
                                GlobalMethods.P_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU5_Dep_textbox.Text);
                            }
                            if (GlobalMethods.landuse[row, col] == 6)
                            {
                                GlobalMethods.infil[row, col] *= System.Convert.ToDouble(landuse_determinator.LU6_Inf_textbox.Text);
                                GlobalMethods.evapotranspiration[row, col] *= System.Convert.ToDouble(landuse_determinator.LU6_Evap_textbox.Text);
                                GlobalMethods.K_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU6_Ero_textbox.Text);
                                GlobalMethods.P_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU6_Dep_textbox.Text);
                            }
                            if (GlobalMethods.landuse[row, col] == 7)
                            {
                                GlobalMethods.infil[row, col] *= System.Convert.ToDouble(landuse_determinator.LU7_Inf_textbox.Text);
                                GlobalMethods.evapotranspiration[row, col] *= System.Convert.ToDouble(landuse_determinator.LU7_Evap_textbox.Text);
                                GlobalMethods.K_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU7_Ero_textbox.Text);
                                GlobalMethods.P_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU7_Dep_textbox.Text);
                            }
                            if (GlobalMethods.landuse[row, col] == 8)
                            {
                                GlobalMethods.infil[row, col] *= System.Convert.ToDouble(landuse_determinator.LU8_Inf_textbox.Text);
                                GlobalMethods.evapotranspiration[row, col] *= System.Convert.ToDouble(landuse_determinator.LU8_Evap_textbox.Text);
                                GlobalMethods.K_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU8_Ero_textbox.Text);
                                GlobalMethods.P_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU8_Dep_textbox.Text);
                            }
                            if (GlobalMethods.landuse[row, col] == 9)
                            {
                                GlobalMethods.infil[row, col] *= System.Convert.ToDouble(landuse_determinator.LU9_Inf_textbox.Text);
                                GlobalMethods.evapotranspiration[row, col] *= System.Convert.ToDouble(landuse_determinator.LU9_Evap_textbox.Text);
                                GlobalMethods.K_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU9_Ero_textbox.Text);
                                GlobalMethods.P_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU9_Dep_textbox.Text);
                            }
                            if (GlobalMethods.landuse[row, col] == 10)
                            {
                                GlobalMethods.infil[row, col] *= System.Convert.ToDouble(landuse_determinator.LU10_Inf_textbox.Text);
                                GlobalMethods.evapotranspiration[row, col] *= System.Convert.ToDouble(landuse_determinator.LU10_Evap_textbox.Text);
                                GlobalMethods.K_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU10_Ero_textbox.Text);
                                GlobalMethods.P_fac[row, col] *= System.Convert.ToDouble(landuse_determinator.LU10_Dep_textbox.Text);
                            }
                        }
                    } //for
                } //for
                  // Debug.WriteLine(" assigned starting values for geomorph  ");
                  // Debug.WriteLine("before initialise soil {0}", GlobalMethods.texture_kg[0, 0, 0, 2]);

                initialise_soil();
                //Debug.WriteLine("after initialise soil {0}", GlobalMethods.texture_kg[0, 0, 0, 2]);
                if (findnegativetexture())
                {
                    Debug.WriteLine("err_ini1");
                    // Debugger.Break(); 
                }

                // displaysoil(0, 0);
                // Debug.WriteLine("Total catchment mass = " + total_catchment_mass());


                //displaysoil(50, 0);
                //writesoil(0, 0);
            } //try
            catch { Debug.WriteLine(" problem assigning starting values to matrices "); }

            if (guiVariables.Fill_sinks_before_checkbox)
            {
                try
                {
                    findsinks();
                    if (GlobalMethods.numberofsinks > 0)
                    {
                        searchdepressions();
                        //define_fillheight_new();
                        for (int row = 0; row < GlobalMethods.nr; row++)
                        {
                            for (int col = 0; col < GlobalMethods.nc; col++)
                            {
                                if (GlobalMethods.dtm[row, col] < GlobalMethods.dtmfill_A[row, col] && GlobalMethods.dtm[row, col] != -9999) { GlobalMethods.dtm[row, col] = GlobalMethods.dtmfill_A[row, col]; }
                            }
                        }
                    }
                }
                catch { Debug.WriteLine(" problem with sink definition "); }
            }

            //Google Earth Preparation
            if (guiVariables.GoogleAnimationCheckbox)
            {
                KML_FILE_NAME = GlobalMethods.Workdir + "\\animation\\animation.kml";
                Debug.WriteLine("creating directory " + GlobalMethods.Workdir + @"\animation");
                Directory.CreateDirectory(GlobalMethods.Workdir + @"\animation");
                startDate = googleBeginDate.Text;
                try { googleTime = System.DateTime.Parse(googleBeginDate.Text); }
                catch { GlobalMethods.input_data_error = true; MessageBox.Show("No valid Google Earth animation start date provided"); }
            }
            // AVI preparation
            if (guiVariables.CheckBoxGenerateAVIFile)
            {
                try
                {
                    aw = new AviWriter();

                    string fname = GlobalMethods.Workdir + @"\" + textBoxAVIFile.Text;
                    Debug.WriteLine("checking for existence of movie file " + fname);
                    //Delete the avi if it already exists, & wait for it to occur
                    if (File.Exists(fname))
                    {
                        File.Delete(fname);

                        if (File.Exists(fname))
                        {
                            System.IO.FileSystemWatcher watcher =
                                new System.IO.FileSystemWatcher(fname);
                            watcher.WaitForChanged(System.IO.WatcherChangeTypes.Deleted);
                        }
                    }
                    bmp = aw.Open(fname, 25, this.Mapwindow.Size.Width,
                        this.Mapwindow.Size.Height);  // <JMW 20041018>
                }
                catch { GlobalMethods.input_data_error = true; MessageBox.Show("Video file can not be created"); }

            }

            // Timeseries preparation
            try
            {
                number_of_outputs = 0;
                if (guiVariables.Timeseries.Timeseries_cell_waterflow_check) { timeseries_order[1] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Timeseries_cell_altitude_check) { timeseries_order[2] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Timeseries_net_ero_check) { timeseries_order[3] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Timeseries_number_dep_check) { timeseries_order[4] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Timeseries_number_erosion_check) { timeseries_order[5] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Timeseries_number_waterflow_check) { timeseries_order[6] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Timeseries_SDR_check) { timeseries_order[7] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Timeseries_total_average_alt_check) { timeseries_order[8] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Timeseries_total_dep_check) { timeseries_order[9] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Timeseries_total_ero_check) { timeseries_order[10] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Timeseries_total_evap_check) { timeseries_order[11] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Timeseries_total_infil_check) { timeseries_order[12] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Timeseries_total_outflow_check) { timeseries_order[13] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Timeseries_total_rain_check) { timeseries_order[14] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Total_phys_weath_checkbox) { timeseries_order[15] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Total_chem_weath_checkbox) { timeseries_order[16] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Total_fine_formed_checkbox) { timeseries_order[17] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Total_fine_eluviated_checkbox) { timeseries_order[18] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Total_mass_bioturbed_checkbox) { timeseries_order[19] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Total_OM_input_checkbox) { timeseries_order[20] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Total_average_soilthickness_checkbox) { timeseries_order[21] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Timeseries_number_soil_thicker_checkbox) { timeseries_order[22] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Timeseries_coarser_checkbox) { timeseries_order[23] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Timeseries_soil_depth_checkbox) { timeseries_order[24] = number_of_outputs; number_of_outputs++; }
                if (guiVariables.Timeseries.Timeseries_soil_mass_checkbox) { timeseries_order[25] = number_of_outputs; number_of_outputs++; }
            }
            catch { Debug.WriteLine("timeseries preparation was unsuccesful"); }

            //if ((Final_output_checkbox.Checked && GlobalMethods.t == guiVariables.End_time) || (Regular_output_checkbox.Checked && (GlobalMethods.t % (int.Parse(guiVariables.Box_years_output)) == 0)))
            //Debug.WriteLine(" successfully ended initialisations  ");
        }

        void initialise_soil_standard()
        {
            double depth_m;
            Debug.WriteLine("initialising soil");
            // At this point, we know the input soildepth at every location (may be zero). 
            // We do not yet know how many layers that corresponds to.
            // If soildepth is not zero, we will calculate the number of layers and assign thicknesses and material to them.
            int soil_layer, texture_class;
            upper_particle_size[0] = Convert.ToDouble(upper_particle_coarse_textbox.Text);
            upper_particle_size[1] = Convert.ToDouble(upper_particle_sand_textbox.Text);
            upper_particle_size[2] = Convert.ToDouble(upper_particle_silt_textbox.Text);
            upper_particle_size[3] = Convert.ToDouble(upper_particle_clay_textbox.Text);
            upper_particle_size[4] = Convert.ToDouble(upper_particle_fine_clay_textbox.Text);
            //calculate bulk density so that we know how much kg of material goes into a layer.  //ART this will go wrong when there are different textures in different locations, but is faster up until that time.
            double coarsefrac = Convert.ToDouble(soildata.coarsebox.Text) / 100;
            double sandfrac = Convert.ToDouble(soildata.sandbox.Text) / 100;
            double siltfrac = Convert.ToDouble(soildata.siltbox.Text) / 100;
            double clayfrac = Convert.ToDouble(soildata.claybox.Text) / 100;
            double fclayfrac = Convert.ToDouble(soildata.fineclaybox.Text) / 100;
            double location_bd;
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    depth_m = 0;
                    if (guiVariables.Creep_testing)
                    {
                        coarsefrac = 0;
                        sandfrac = 1;
                        siltfrac = 0;
                        clayfrac = 0;
                        fclayfrac = 0;

                    }

                    if (GlobalMethods.soildepth_m[row, col] == 0)
                    {

                        for (soil_layer = 0; soil_layer < GlobalMethods.max_soil_layers; soil_layer++)
                        {
                            for (texture_class = 0; texture_class < GlobalMethods.n_texture_classes; texture_class++)
                            {
                                GlobalMethods.texture_kg[row, col, soil_layer, texture_class] = 0;
                            }
                            GlobalMethods.young_SOM_kg[row, col, soil_layer] = 0;
                            GlobalMethods.old_SOM_kg[row, col, soil_layer] = 0;
                            GlobalMethods.bulkdensity[row, col, soil_layer] = 0;
                            GlobalMethods.layerthickness_m[row, col, soil_layer] = -1;
                        }
                    }
                    else
                    {

                        //now assign thicknesses and material to layer.
                        double available_soildepth = GlobalMethods.soildepth_m[row, col];
                        soil_layer = 0;
                        while (available_soildepth > 0)
                        {
                            // 0-50 cm    min 2.5   insteek 5    maximum 10 cm       n=10    bovenste laag geen minimum (sediment HOEFT niet meteen weggemiddeld te worden - pas als nodig)
                            // 50-200 cm  min 10    insteek 15    maximum 50 cm      n=10
                            // daarna     min 50    insteek 100  geen max            n=5
                            // If GlobalMethods.max_soil_layers is smaller than the sum of the perfect layers in each of the three ' packages' , then we simply make the lowest layer very thick.
                            //if (soil_layer < 10 && soil_layer < GlobalMethods.max_soil_layers - 1)
                            /*
                            if (soil_layer < 40 && soil_layer < GlobalMethods.max_soil_layers - 1)
                            {
                                if (available_soildepth > 0.05)
                                {
                                    GlobalMethods.layerthickness_m[row, col, soil_layer] = 0.05;
                                    available_soildepth -= 0.05;
                                }
                                else
                                {
                                    GlobalMethods.layerthickness_m[row, col, soil_layer] = available_soildepth;
                                    available_soildepth = 0;
                                }
                            }
                            */
                            if (soil_layer < 10 && soil_layer < GlobalMethods.max_soil_layers - 1)
                            {
                                if (available_soildepth > 0.05) //
                                {
                                    GlobalMethods.layerthickness_m[row, col, soil_layer] = 0.05; // 
                                    available_soildepth -= 0.05;
                                }
                                else
                                {
                                    GlobalMethods.layerthickness_m[row, col, soil_layer] = available_soildepth;
                                    available_soildepth = 0;
                                }
                            }
                            if (soil_layer > 9 && soil_layer < 20 && soil_layer < GlobalMethods.max_soil_layers - 1)
                            {
                                if (available_soildepth > 0.15) // was 0.25
                                {
                                    GlobalMethods.layerthickness_m[row, col, soil_layer] = 0.15; // was 0.15
                                    available_soildepth -= 0.15;
                                }
                                else
                                {
                                    GlobalMethods.layerthickness_m[row, col, soil_layer] = available_soildepth;
                                    available_soildepth = 0;
                                }
                            }
                            if (soil_layer > 19 && soil_layer < GlobalMethods.max_soil_layers && soil_layer < GlobalMethods.max_soil_layers - 1) // Rest
                            {
                                if (available_soildepth > 0.5) // was 1
                                {
                                    GlobalMethods.layerthickness_m[row, col, soil_layer] = 0.5; // was 1
                                    available_soildepth -= 0.5; // was 1
                                }
                                else
                                {
                                    GlobalMethods.layerthickness_m[row, col, soil_layer] = available_soildepth;
                                    available_soildepth = 0;
                                }
                            }

                            if (soil_layer == GlobalMethods.max_soil_layers - 1)
                            {
                                GlobalMethods.layerthickness_m[row, col, soil_layer] = available_soildepth;
                                available_soildepth = 0;
                            }

                            if (GlobalMethods.layerthickness_m[row, col, soil_layer] != 0)
                            {
                                depth_m += GlobalMethods.layerthickness_m[row, col, soil_layer] / 2;
                                location_bd = bulk_density_calc(coarsefrac, sandfrac, siltfrac, clayfrac, fclayfrac, 0, 0, depth_m);
                                depth_m += GlobalMethods.layerthickness_m[row, col, soil_layer] / 2;
                                GlobalMethods.texture_kg[row, col, soil_layer, 0] = location_bd * GlobalMethods.layerthickness_m[row, col, soil_layer] * coarsefrac * GlobalMethods.dx * GlobalMethods.dx;   //  kg = kg/m3 * m * kg/kg * m * m
                                GlobalMethods.texture_kg[row, col, soil_layer, 1] = location_bd * GlobalMethods.layerthickness_m[row, col, soil_layer] * sandfrac * GlobalMethods.dx * GlobalMethods.dx;
                                GlobalMethods.texture_kg[row, col, soil_layer, 2] = location_bd * GlobalMethods.layerthickness_m[row, col, soil_layer] * siltfrac * GlobalMethods.dx * GlobalMethods.dx;
                                GlobalMethods.texture_kg[row, col, soil_layer, 3] = location_bd * GlobalMethods.layerthickness_m[row, col, soil_layer] * clayfrac * GlobalMethods.dx * GlobalMethods.dx;
                                GlobalMethods.texture_kg[row, col, soil_layer, 4] = location_bd * GlobalMethods.layerthickness_m[row, col, soil_layer] * fclayfrac * GlobalMethods.dx * GlobalMethods.dx;
                                GlobalMethods.bulkdensity[row, col, soil_layer] = location_bd;

                                if (guiVariables.Decalcification_checkbox)
                                {
                                    CO3_kg[row, col, soil_layer] = (location_bd * GlobalMethods.layerthickness_m[row, col, soil_layer] * GlobalMethods.dx * GlobalMethods.dx) * Convert.ToDouble(ini_CaCO3_content.Text) * 40.08 / (40.08 + 60.01); // calculate total CO3: total mass * fraction of soil * fraction of CaCO3 molecule
                                }

                            }
                            if (guiVariables.Creep_testing)
                            {
                                sandfrac -= 0.05;
                                clayfrac += 0.05;

                                sandfrac = Math.Max(sandfrac, 0);
                                clayfrac = Math.Min(clayfrac, 1);
                            }

                            soil_layer++;

                        } // end availabke soil depth > 0
                    } // end else 

                } // end col
            } // end row
              // Debug.WriteLine("initialised soil");

        }

        void initialise_soil()
        {
            double depth_m;
            //Debug.WriteLine("initialising soil");
            // At this point, we know the input soildepth at every location (may be zero). 
            // We do not yet know how many layers that corresponds to.
            // If soildepth is not zero, we will calculate the number of layers and assign thicknesses and material to them.
            int soil_layer, texture_class;
            upper_particle_size[0] = Convert.ToDouble(upper_particle_coarse_textbox.Text);
            upper_particle_size[1] = Convert.ToDouble(upper_particle_sand_textbox.Text);
            upper_particle_size[2] = Convert.ToDouble(upper_particle_silt_textbox.Text);
            upper_particle_size[3] = Convert.ToDouble(upper_particle_clay_textbox.Text);
            upper_particle_size[4] = Convert.ToDouble(upper_particle_fine_clay_textbox.Text);
            //calculate bulk density so that we know how much kg of material goes into a layer.  //ART this will go wrong when there are different textures in different locations, but is faster up until that time.
            double coarsefrac = Convert.ToDouble(soildata.coarsebox.Text) / 100;
            double sandfrac = Convert.ToDouble(soildata.sandbox.Text) / 100;
            double siltfrac = Convert.ToDouble(soildata.siltbox.Text) / 100;
            double clayfrac = Convert.ToDouble(soildata.claybox.Text) / 100;
            double fclayfrac = Convert.ToDouble(soildata.fineclaybox.Text) / 100;
            double location_bd;
            double dz_standard = 0.1;
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    depth_m = 0;
                    if (guiVariables.Creep_testing)
                    {
                        coarsefrac = 0;
                        sandfrac = 1;
                        siltfrac = 0;
                        clayfrac = 0;
                        fclayfrac = 0;

                    }

                    if (GlobalMethods.soildepth_m[row, col] == 0)
                    {

                        for (soil_layer = 0; soil_layer < GlobalMethods.max_soil_layers; soil_layer++)
                        {
                            for (texture_class = 0; texture_class < GlobalMethods.n_texture_classes; texture_class++)
                            {
                                GlobalMethods.texture_kg[row, col, soil_layer, texture_class] = 0;
                            }
                            GlobalMethods.young_SOM_kg[row, col, soil_layer] = 0;
                            GlobalMethods.old_SOM_kg[row, col, soil_layer] = 0;
                            GlobalMethods.bulkdensity[row, col, soil_layer] = 0;
                            GlobalMethods.layerthickness_m[row, col, soil_layer] = -1;
                        }
                    }
                    else
                    {

                        //now assign thicknesses and material to layer.
                        double available_soildepth = GlobalMethods.soildepth_m[row, col];
                        soil_layer = 0;
                        while (available_soildepth > 0)
                        {
                            if (available_soildepth > dz_standard)
                            {
                                GlobalMethods.layerthickness_m[row, col, soil_layer] = dz_standard;
                                available_soildepth -= dz_standard;
                            }
                            else
                            {
                                GlobalMethods.layerthickness_m[row, col, soil_layer] = available_soildepth;
                                available_soildepth = 0;
                            }

                            if (soil_layer == GlobalMethods.max_soil_layers - 1)
                            {
                                GlobalMethods.layerthickness_m[row, col, soil_layer] = available_soildepth;
                                available_soildepth = 0;
                            }

                            if (GlobalMethods.layerthickness_m[row, col, soil_layer] != 0)
                            {
                                depth_m += GlobalMethods.layerthickness_m[row, col, soil_layer] / 2;
                                location_bd = bulk_density_calc(coarsefrac, sandfrac, siltfrac, clayfrac, fclayfrac, 0, 0, depth_m);
                                depth_m += GlobalMethods.layerthickness_m[row, col, soil_layer] / 2;
                                GlobalMethods.texture_kg[row, col, soil_layer, 0] = location_bd * GlobalMethods.layerthickness_m[row, col, soil_layer] * coarsefrac * GlobalMethods.dx * GlobalMethods.dx;   //  kg = kg/m3 * m * kg/kg * m * m
                                GlobalMethods.texture_kg[row, col, soil_layer, 1] = location_bd * GlobalMethods.layerthickness_m[row, col, soil_layer] * sandfrac * GlobalMethods.dx * GlobalMethods.dx;
                                GlobalMethods.texture_kg[row, col, soil_layer, 2] = location_bd * GlobalMethods.layerthickness_m[row, col, soil_layer] * siltfrac * GlobalMethods.dx * GlobalMethods.dx;
                                GlobalMethods.texture_kg[row, col, soil_layer, 3] = location_bd * GlobalMethods.layerthickness_m[row, col, soil_layer] * clayfrac * GlobalMethods.dx * GlobalMethods.dx;
                                GlobalMethods.texture_kg[row, col, soil_layer, 4] = location_bd * GlobalMethods.layerthickness_m[row, col, soil_layer] * fclayfrac * GlobalMethods.dx * GlobalMethods.dx;
                                GlobalMethods.bulkdensity[row, col, soil_layer] = location_bd;

                                if (guiVariables.Decalcification_checkbox)
                                {
                                    CO3_kg[row, col, soil_layer] = (location_bd * GlobalMethods.layerthickness_m[row, col, soil_layer] * GlobalMethods.dx * GlobalMethods.dx) * Convert.ToDouble(ini_CaCO3_content.Text) * 40.08 / (40.08 + 60.01); // calculate total CO3: total mass * fraction of soil * fraction of CaCO3 molecule
                                }

                            }
                            if (guiVariables.Creep_testing)
                            {
                                sandfrac -= 0.05;
                                clayfrac += 0.05;

                                sandfrac = Math.Max(sandfrac, 0);
                                clayfrac = Math.Min(clayfrac, 1);
                            }

                            soil_layer++;

                        } // end availabke soil depth > 0
                    } // end else 
                    GlobalMethods.soildepth_m[row, col] = total_soil_thickness(row, col);

                } // end col
            } // end row
            //Debug.WriteLine("initialised soil");

        } // anngepast voor standaard diktes


        void initialise_every_till()
        {
            if (guiVariables.Check_time_till_fields && !guiVariables.Check_space_till_fields)
            {
                for (int row = 0; row < GlobalMethods.nr; row++)
                {
                    for (int col = 0; col < GlobalMethods.nc; col++)
                    {
                        GlobalMethods.tillfields[row, col] = 1 * till_record[GlobalMethods.t];
                    }
                }

            }
        }

        void initialise_every()                         //fills the inputgrids with values
        {
            int corrected_t;
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    // time runs from 1 to guiVariables.End_time - compensate for that when taking values from records
                    // also compensate for records shorter than guiVariables.End_time
                    if (check_time_rain.Checked)
                    {
                        corrected_t = GlobalMethods.t;
                        while (corrected_t > rainfall_record.Length) { corrected_t -= rainfall_record.Length; }

                        rain_value_m = 0.001 * rainfall_record[corrected_t]; //from mm (in record) to m (LORICA)   // mvdm -1 weggehaald van corrected_t, leidde tot OutOfRange errors
                                                                             // changed GlobalMethods.rain[row, col] to rain_value_m, due to errors, this is not spatial, but temporal variation
                                                                             //this should be improved for when rainfall is not also spatially variable //ArT
                    }
                    if (check_time_infil.Checked)
                    {
                        corrected_t = GlobalMethods.t;
                        while (corrected_t > infil_record.Length) { corrected_t -= infil_record.Length; }
                        infil_value_m = 0.001 * infil_record[corrected_t];
                    }
                    if (check_time_evap.Checked)
                    {
                        corrected_t = GlobalMethods.t;
                        while (corrected_t > evap_record.Length) { corrected_t -= evap_record.Length; }
                        evap_value_m = 0.001 * evap_record[corrected_t];
                    }

                    if (check_time_T.Checked)
                    {
                        corrected_t = GlobalMethods.t;
                        while (corrected_t > temp_record.Length) { corrected_t -= temp_record.Length; }

                        rain_value_m = 0.001 * temp_record[corrected_t]; //from mm (in record) to m (LORICA)   // mvdm -1 weggehaald van corrected_t, leidde tot OutOfRange errors
                        temp_value_C = temp_record[corrected_t];
                        // changed GlobalMethods.rain[row, col] to rain_value_m, due to errors, this is not spatial, but temporal variation
                        //this should be improved for when rainfall is not also spatially variable //ArT
                    }




                    if (annual_output_checkbox.Checked)
                    {
                        GlobalMethods.dz_soil[row, col] = 0;
                        if (Creep_Checkbox.Checked) { GlobalMethods.sum_creep_grid[row, col] = 0; GlobalMethods.creep[row, col] = 0; }
                        if (treefall_checkbox.Checked) { GlobalMethods.dz_treefall[row, col] = 0; GlobalMethods.treefall_count[row, col] = 0; }
                        if (Solifluction_checkbox.Checked) { GlobalMethods.sum_solifluction[row, col] = 0; }
                        if (Water_ero_checkbox.Checked) { GlobalMethods.sum_water_erosion[row, col] = 0; }
                        if (Biological_weathering_checkbox.Checked) { GlobalMethods.sum_biological_weathering[row, col] = 0; }
                        if (Frost_weathering_checkbox.Checked) { GlobalMethods.sum_frost_weathering[row, col] = 0; }
                        if (Tillage_checkbox.Checked) { GlobalMethods.sum_tillage[row, col] = 0; total_sum_tillage = 0; }
                        if (GlobalMethods.soildepth_m[row, col] < 0.0) { soildepth_error += GlobalMethods.soildepth_m[row, col]; GlobalMethods.soildepth_m[row, col] = 0; }
                    }
                    if (GlobalMethods.soildepth_m[row, col] < 0.0) { soildepth_error += GlobalMethods.soildepth_m[row, col]; GlobalMethods.soildepth_m[row, col] = 0; }
                    if (Water_ero_checkbox.Checked) { GlobalMethods.waterflow_m3[row, col] = 0.0; }

                } //for
            } //for

            if (fill_sinks_during_checkbox.Checked)
            {
                findsinks();
                if (GlobalMethods.numberofsinks > 0)
                {
                    searchdepressions();
                    //define_fillheight_new();
                    for (int row = 0; row < GlobalMethods.nr; row++)
                    {
                        for (int col = 0; col < GlobalMethods.nc; col++)
                        {
                            if (GlobalMethods.dtm[row, col] < GlobalMethods.dtmfill_A[row, col] && GlobalMethods.dtm[row, col] != -9999) { GlobalMethods.dtm[row, col] = GlobalMethods.dtmfill_A[row, col]; }
                        }
                    }
                }
            }
            //Debug.WriteLine("initialised every");
        }

        #endregion

        #region calibration code

        void makerecords(string filename)
        {
            string FILE_NAME = filename;
            string input;
            if (!File.Exists(FILE_NAME))
            {
                MessageBox.Show("No such data file " + FILE_NAME);
                GlobalMethods.input_data_error = true;
                return;
            }
            Debug.WriteLine("reading " + filename + " into record ");
            StreamReader sr = File.OpenText(FILE_NAME);

            //read first line: number of timesteps
            input = sr.ReadLine();
            int recordsize = 0;
            try { recordsize = System.Convert.ToInt32(input); }
            catch
            {
                MessageBox.Show("Wrong value " + input + " in first line of record " + FILE_NAME);
                GlobalMethods.input_data_error = true;
                return;
            }
            if (check_time_rain.Checked) { rainfall_record = new int[recordsize]; }
            if (check_time_evap.Checked) { evap_record = new int[recordsize]; }
            if (check_time_infil.Checked) { infil_record = new int[recordsize]; }
            if (check_time_T.Checked) { temp_record = new int[recordsize]; }
            if (guivariables.Check_time_till_fields) { till_record = new int[recordsize]; }

            memory_records = true;
        }

        void makedailyrecords(string filename)
        {
            string FILE_NAME = filename;
            string input;
            if (!File.Exists(FILE_NAME))
            {
                MessageBox.Show("No such data file " + FILE_NAME);
                GlobalMethods.input_data_error = true;
                return;
            }
            Debug.WriteLine("reading " + filename + " into record ");
            StreamReader sr = File.OpenText(FILE_NAME);

            //read first line: number of timesteps
            input = sr.ReadLine();
            int recordsize = 0;
            try { recordsize = System.Convert.ToInt32(input); }
            catch
            {
                MessageBox.Show("Wrong value " + input + " in first line of record " + FILE_NAME);
                GlobalMethods.input_data_error = true;
                return;
            }


            guiVariables.P_all = new int[recordsize];
            guiVariables.ET0_all = new int[recordsize];
            guiVariables.D_all = new int[recordsize];
            guiVariables.Tavg_all = new int[recordsize];
            guiVariables.Tmin_all = new int[recordsize];
            guiVariables.Tmax_all = new int[recordsize];


            memory_records_d = true;
        }

        private void calib_calculate_maxruns(int calibparacount)
        {
            //this code calculates the total number of runs needed when calibrating
            string calibration_ratio_string = calibration_ratios_textbox.Text;
            string[] ratiowords = calibration_ratio_string.Split(';');
            int ratio;
            for (ratio = 0; ratio < ratiowords.Length; ratio++)
                for (int par = 0; par < calibparacount; par++)
                {
                    try
                    {
                        calib_ratios[par, ratio] = Convert.ToDouble(ratiowords[ratio]);
                    }
                    catch { GlobalMethods.input_data_error = true; MessageBox.Show("Calibration ratio input error"); }
                }
            try { GUcalib_levels = Convert.ToInt32(calibration_levels_textbox.Text); }
            catch { GlobalMethods.input_data_error = true; MessageBox.Show("Calibration iterations must be an integer"); }
            maxruns = calib_levels * Convert.ToInt32(Math.Pow(ratiowords.Length, calibparacount));
            Debug.WriteLine(" the number of runs for calibration will be " + maxruns);
        }

        private void calib_shift_and_zoom(int para_number, double zoom_factor, double orig_par_value)
        {
            //this code iinds out whether the best parameter value was on the edge or inside the range explored. Then shifts and zooms out or in , depending
            try
            {
                Debug.WriteLine(" para number " + para_number);
                Debug.WriteLine(" best_parameter value " + best_parameters[para_number]);
                double mid_ratio = 0;
                if (calib_ratios.GetLength(1) % 2 == 0) { mid_ratio = (calib_ratios[para_number, Convert.ToInt32(calib_ratios.GetLength(1) / 2) - 1] + calib_ratios[para_number, Convert.ToInt32((calib_ratios.GetLength(1) / 2))]) / 2; }
                else { mid_ratio = calib_ratios[para_number, Convert.ToInt32(calib_ratios.GetLength(1) / 2 - 0.5)]; }
                Debug.WriteLine("mid ratio is " + mid_ratio);
                Double best_ratio = best_parameters[para_number] / orig_par_value;
                Debug.WriteLine("best ratio is " + best_ratio);
                if (best_parameters[para_number] == calib_ratios[para_number, 0] * orig_par_value | best_parameters[para_number] == calib_ratios[para_number, calib_ratios.GetLength(1) - 1] * orig_par_value)
                {
                    //the best parameter ratio (and thus value) was on the edge of the range. We must shift our range sideways (we keep the same ratio between upper and lower ratio - are you still with me?)
                    Debug.WriteLine(" currentpara value was on edge of range");

                    for (int ratio = 0; ratio < calib_ratios.GetLength(1); ratio++)
                    {
                        Debug.WriteLine(" setting ratio " + calib_ratios[para_number, ratio] + " to " + calib_ratios[para_number, ratio] * (best_ratio / mid_ratio));
                        calib_ratios[para_number, ratio] = calib_ratios[para_number, ratio] * (best_ratio / mid_ratio);
                    }
                }
                else
                {
                    //the best parameter ratio (and thus value) NOT on the edge of the range. We must shift to the best observed value and then zoom IN
                    Debug.WriteLine(" currentpara value was NOT on edge of range");
                    for (int ratio = 0; ratio < calib_ratios.GetLength(1); ratio++)
                    {
                        Debug.WriteLine(" setting ratio " + calib_ratios[para_number, ratio] + " to " + ((best_ratio - calib_ratios[para_number, ratio]) / zoom_factor));
                        calib_ratios[para_number, ratio] += (best_ratio - calib_ratios[para_number, ratio]) / zoom_factor;
                    }
                }
            }
            catch { Debug.WriteLine(" problem adapting parameters and ratios "); }
        }

        private void calib_prepare_report()
        {
            //this code prepares a calibration report
            //it opens and writes headers for a text file on disk
            string FILENAME = GlobalMethods.Workdir + "\\calibration.log";
            using (StreamWriter sw = new StreamWriter(FILENAME))
            {
                try
                {
                    sw.Write("run objective_function_value ");
                    //USER INPUT NEEDED IN FOLLOWING LINE: ENTER THE CALIBRATION PARAMETER NAMES 
                    //THEY WILL BE HEADERS IN THE CALIBRATION REPORT
                    sw.WriteLine(" erodibility_K conv_fac");
                }
                catch { Debug.WriteLine(" issue with writing the header of the calibration log file"); }
            }
            Debug.WriteLine(" calib tst - calib_prepare_rep - added first line to file");
        }

        private void calib_update_report(double objective_fnct_result)
        {
            //this code updates a calibration report
            //it writes parameters and objective function outcomes to disk
            string FILENAME = GlobalMethods.Workdir + "\\calibration.log";
            using (StreamWriter sw = File.AppendText(FILENAME))
            {
                try
                {
                    //USER INPUT NEEDED IN FOLLOWING LINE: ENTER THE CALIBRATION PARAMETERS 
                    sw.WriteLine(run_number + " " + objective_fnct_result + " " + advection_erodibility + " " + conv_fac);
                }
                catch { Debug.WriteLine(" issue with writing a line in the calibration log file"); }
            }
            Debug.WriteLine(" calib tst - calib_update_rep - added line to file");
        }

        private void calib_finish_report()
        {
            //this code closes a calibration report
            //it writes the parameters for the best run to disk
            //CALIB_USER : Change the number of parameters referenced (now two)
            try
            {
                string FILENAME = GlobalMethods.Workdir + "\\calibration.log";
                using (StreamWriter sw = File.AppendText(FILENAME))
                {
                    sw.WriteLine(best_run + " " + best_error + " " + best_parameters[0] + " " + best_parameters[1]);
                    Debug.WriteLine(" best run was " + best_run + " with error " + best_error + "m3");
                }
                Debug.WriteLine(" calib tst - calib_finish_rep - wrote final line and closed file");
            }
            catch
            {
                Debug.WriteLine(" calib tst - calib_finish_rep - FAILED to write file");
            }
        }

        private double calib_objective_function()
        {
            //this code calculates the value of the objective function during calibration and is user-specified. 
            //calibration looks to minimize the value of the objective function by varying parameter values
            //CALIB_USER
            //example for Luxembourg: we want to simulate the correct amount of erosion, over the entire slope
            //Xia, number needs to be adapted
            double simulated_ero_m3 = 0;
            double simulated_ero_kg_m2_y = 0;
            double known_ero_kg_m2_y = 0.0313;
            double total_bulk_density = 0;
            double average_bulk_density = 0;
            int objective_function_cells = 0;
            for (int row = 0; row < GlobalMethods.nr; row++)
            {
                for (int col = 0; col < GlobalMethods.nc; col++)
                {
                    if (GlobalMethods.dtm[row, col] != -9999)
                    {
                        simulated_ero_m3 -= GlobalMethods.sum_water_erosion[row, col] * GlobalMethods.dx * GlobalMethods.dx;
                        total_bulk_density += GlobalMethods.bulkdensity[row, col, 0];
                        objective_function_cells++;
                    }
                }
            }
            average_bulk_density = total_bulk_density / objective_function_cells;
            simulated_ero_kg_m2_y = average_bulk_density * simulated_ero_m3 / guiVariables.End_time / (objective_function_cells * GlobalMethods.dx * GlobalMethods.dx);
            ;
            Debug.WriteLine(" calib tst - calib_objective_function - error is " + Math.Abs(known_ero_kg_m2_y - simulated_ero_kg_m2_y) + "kg per m2 per year");
            return Math.Abs(known_ero_kg_m2_y - simulated_ero_kg_m2_y);
        }

        private void calib_update_best_paras()
        {
            //this code updates the recorded set of parameter values that gives the best score for the objective function
            //USERS have to update code here to reflect the parameters they actually vary
            Debug.WriteLine(" updating parameter set for best scored run");
            // add/change lines below
            best_parameters[0] = advection_erodibility;
            //best_parameters[1] = conv_fac;
            Debug.WriteLine(" best erodib " + best_parameters[0]);
            //Debug.WriteLine(" best conv_fac " + best_parameters[1]);
        }

        #endregion
    }
}
