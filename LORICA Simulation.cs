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

namespace LORICA4
{
    class Simulation
    {
        Stopwatch stopwatch;
        GUIVariables guiVariables;
        public Simulation(GUIVariables gv)
        {
            guiVariables = gv;
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
                workdir = separate[0];
                Debug.WriteLine("storing results in " + workdir);
                System.IO.Directory.CreateDirectory(workdir);
                input_data_error = false;

                try { end_time = int.Parse(guiVariables.Number_runs_textbox); }
                catch { input_data_error = true; MessageBox.Show("Invalid number of years"); }
                try { ntr = System.Convert.ToInt32(end_time); }     // WVG initialise ntr: number of rows in timeseries matrix   
                catch (OverflowException)
                {
                    MessageBox.Show("number of timesteps is outside the range of the Int32 type.");
                }
                //WVG initialise ntr, nr of timesteps, can be changed to nr of output timesteps
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
                if (guiVariables.soil_carbon_cycle_checkbox) //:)
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
                    catch { input_data_error = true; MessageBox.Show("value for google save interval is not valid"); }

                    if (guiVariables.UTMgridcheckbox)
                    {
                        try { test = System.Convert.ToInt32(guiVariables.UTMzonebox); }
                        catch { input_data_error = true; MessageBox.Show("value for UTM zone is not valid"); }
                    }

                    if (end_time < save_interval2 && guiVariables.GoogleAnimationCheckbox) { input_data_error = true; MessageBox.Show("value for google save interval cannot be larger than number of runs "); }
                    if (end_time < int.Parse(guiVariables.Saveintervalbox) && guiVariables.CheckBoxGenerateAVIFile) { input_data_error = true; MessageBox.Show("value for video interval cannot be larger than number of runs "); }


                    //WATER EROSION AND DEPOSITION PARAMETERS
                    if (water_ero_active)
                    {
                        try { m = double.Parse(guiVariables.Parameter_m_textbox); }
                        catch { input_data_error = true; MessageBox.Show("value for parameter m is not valid"); }                      // Kirkby's m and n factors for increasing
                        try { n = double.Parse(guiVariables.Parameter_n_textbox); }
                        catch { input_data_error = true; MessageBox.Show("value for parameter n is not valid"); }                   // sheet, wash, overland, gully to river flow
                        try { conv_fac = double.Parse(guiVariables.Parameter_m_textbox); }
                        catch { input_data_error = true; MessageBox.Show("value for parameter p is not valid"); }
                        try { advection_erodibility = double.Parse(guiVariables.Parameter_K_textbox); }
                        catch { input_data_error = true; MessageBox.Show("value for parameter K is not valid"); }
                        try { bio_protection_constant = double.Parse(guiVariables.Bio_protection_constant_textbox); }
                        catch { input_data_error = true; MessageBox.Show("value for parameter P is not valid"); }
                        try { rock_protection_constant = double.Parse(guiVariables.Rock_protection_constant_textbox); }
                        catch { input_data_error = true; MessageBox.Show("value for parameter P is not valid"); }
                        try { constant_selective_transcap = double.Parse(guiVariables.Selectivity_constant_textbox); }
                        catch { input_data_error = true; MessageBox.Show("value for parameter P is not valid"); }
                        try { erosion_threshold_kg = double.Parse(guiVariables.Erosion_threshold_textbox); }
                        catch { input_data_error = true; MessageBox.Show("value for parameter P is not valid"); }
                    }

                    //TILLAGE PARAMETERS
                    if (tillage_active)
                    {
                        try { plough_depth = double.Parse(guiVariables.Parameter_ploughing_depth_textbox); }
                        catch { input_data_error = true; MessageBox.Show("value for parameter plough depth is not valid"); }
                        try { tilc = double.Parse(guiVariables.Parameter_tillage_constant_textbox); }
                        catch { input_data_error = true; MessageBox.Show("value for parameter tillage constant is not valid"); }
                    }

                    //CREEP PARAMETER
                    if (creep_active)
                    {
                        try { conv_fac = double.Parse(guiVariables.Parameter_m_textbox); }
                        catch { input_data_error = true; MessageBox.Show("value for parameter p is not valid"); }
                        try { diffusivity_creep = double.Parse(guiVariables.Parameter_diffusivity_textbox); }
                        catch { input_data_error = true; MessageBox.Show("value for parameter diffusivity is not valid"); }
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
                        catch { input_data_error = true; MessageBox.Show("value for parameter P0 is not valid"); }
                        try { k1 = double.Parse(guiVariables.Parameter_k1_textbox); }
                        catch { input_data_error = true; MessageBox.Show("value for parameter k1 is not valid"); }
                        try { k2 = double.Parse(guiVariables.Parameter_k2_textbox); }
                        catch { input_data_error = true; MessageBox.Show("value for parameter k2 is not valid"); }
                        try { Pa = double.Parse(guiVariables.Parameter_Pa_textbox); }
                        catch { input_data_error = true; MessageBox.Show("value for parameter Pa is not valid"); }
                    }

                    //Tilting parameters
                    if (tilting_active)
                    {
                        if (guiVariables.Radio_tilt_col_zero) { tilt_location = 0; }
                        if (guiVariables.Radio_tilt_row_zero) { tilt_location = 1; }
                        if (guiVariables.Radio_tilt_col_max) { tilt_location = 2; }
                        if (guiVariables.Radio_tilt_row_max) { tilt_location = 3; }
                        try { tilt_intensity = double.Parse(guiVariables.Tilting_rate_textbox); }
                        catch { input_data_error = true; MessageBox.Show("value for parameter tilting rate is not valid"); }
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
                            catch { input_data_error = true; MessageBox.Show("value for parameter tilting rate is not valid"); }
                        }
                        if (lift_type == 1)
                        {
                            try { lift_location = int.Parse(guiVariables.Text_lift_row_more); }
                            catch { input_data_error = true; MessageBox.Show("value for parameter tilting rate is not valid"); }
                        }
                        if (lift_type == 2)
                        {
                            try { lift_location = int.Parse(guiVariables.Text_lift_col_less); }
                            catch { input_data_error = true; MessageBox.Show("value for parameter tilting rate is not valid"); }
                        }
                        if (lift_type == 3)
                        {
                            try { lift_location = int.Parse(guiVariables.Text_lift_col_more); }
                            catch { input_data_error = true; MessageBox.Show("value for parameter tilting rate is not valid"); }
                        }
                        try { lift_intensity = double.Parse(guiVariables.Uplift_rate_textbox); }
                        catch { input_data_error = true; MessageBox.Show("value for parameter tilting rate is not valid"); }
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
                            input_data_error = true; Debug.WriteLine("problem reading parameters for pysical weathering");
                        }
                    }


                    //SOIL CHEMICAL WEATHERING PARAMETERS
                    if (soil_chem_weath_active)
                    {
                        try
                        {
                            chemical_weathering_constant = Convert.ToDouble(chem_weath_rate_constant_textbox.Text);
                            Cthree = Convert.ToDouble(chem_weath_depth_constant_textbox.Text);
                            Cfour = Convert.ToDouble(chem_weath_specific_coefficient_textbox.Text);
                            specific_area[0] = Convert.ToDouble(specific_area_coarse_textbox.Text);
                            specific_area[1] = Convert.ToDouble(specific_area_sand_textbox.Text);
                            specific_area[2] = Convert.ToDouble(specific_area_silt_textbox.Text);
                            specific_area[3] = Convert.ToDouble(specific_area_clay_textbox.Text);
                            specific_area[4] = Convert.ToDouble(specific_area_fine_clay_textbox.Text);
                            neoform_constant = Convert.ToDouble(clay_neoform_constant_textbox.Text);
                            Cfive = Convert.ToDouble(clay_neoform_C1_textbox.Text);
                            Csix = Convert.ToDouble(clay_neoform_C2_textbox.Text);
                            Debug.WriteLine("succesfully read parameters for chemical weathering");
                        }
                        catch
                        {
                            input_data_error = true; Debug.WriteLine("problem reading parameters for chemical weathering");
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
                            input_data_error = true; Debug.WriteLine("problem reading parameters for clay dynamics");
                        }
                        if (guiVariables.CT_depth_decay_checkbox)
                        {
                            try
                            {
                                ct_depthdec = Convert.ToDouble(guiVariables.Ct_depth_decay);
                            }
                            catch
                            {
                                input_data_error = true; Debug.WriteLine("problem reading depth decay parameter for clay dynamics");
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
                            input_data_error = true; Debug.WriteLine("problem reading parameters for bioturbation");
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
                            input_data_error = true; Debug.WriteLine("problem reading parameters for carbon cycle");
                        }
                    }

                    try
                    {
                        filename = dtmfilename;             //for directory input
                        dtm_file(filename);                 // from dtm_file(), almost all memory for the model is claimed
                    }
                    catch { Debug.WriteLine(" failed to initialise dtm "); }

                    //LARGEST THING IN HERE. RUNS FOREVER (very long time)
                    if (input_data_error == false)
                    {
                        try
                        {

                            //Debug.WriteLine("reading general values");
                            if (check_space_soildepth.Checked != true)
                            {
                                try { soildepth_value = double.Parse(soildepth_constant_value_box.Text); }
                                catch { MessageBox.Show("value for parameter soildepth is not valid"); }
                            }
                            if (check_space_landuse.Checked != true && check_time_landuse.Checked != true)
                            {
                                try { landuse_value = int.Parse(landuse_constant_value_box.Text); }
                                catch { MessageBox.Show("value for parameter landuse is not valid"); }
                            }
                            if (check_space_evap.Checked != true && check_time_evap.Checked != true)
                            {
                                try { evap_value_m = double.Parse(evap_constant_value_box.Text); }
                                catch { MessageBox.Show("value for parameter evapotranspiration is not valid"); }
                            }
                            if (check_space_infil.Checked != true && check_time_infil.Checked != true)
                            {
                                try { infil_value_m = double.Parse(infil_constant_value_box.Text); }
                                catch { MessageBox.Show("value for parameter infiltration is not valid"); }
                            }

                            if (check_space_rain.Checked != true && check_time_rain.Checked != true)
                            {
                                try { rain_value_m = double.Parse(rainfall_constant_value_box.Text); }
                                catch { MessageBox.Show("value for parameter rainfall is not valid"); }
                            }

                            if (check_time_T.Checked != true)
                            {
                                try { temp_value_C = int.Parse(temp_constant_value_box.Text); }
                                catch { MessageBox.Show("value for parameter temperature is not valid"); }
                            }

                        }
                        catch { MessageBox.Show("there was a problem reading input values"); input_data_error = true; }
                        // Debug.WriteLine("initialising non-general inputs");
                        try { initialise_once(); } // reading input files
                        catch { MessageBox.Show("there was a problem reading input files "); input_data_error = true; }


                        //CALIB_USER: multiply parameter values with current ratio
                        //Note the correspondence between the formulas. Change only 1 value for additional parameters!
                        if (Calibration_button.Checked == true)
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

                        timeseries_matrix = new double[System.Convert.ToInt32(end_time), number_of_outputs];
                        if (input_data_error == false)
                        {
                            try
                            {
                                /*tabControl1.Visible = false;  
                                Mapselector.Enabled = true; 
                                try {View_tabs_checkbox.Checked = false;}
                                catch { Debug.WriteLine(" failed to set view_tabs to unchecked "); }
                                graphicToGoogleEarthButton.Visible = true;
                                graphics_scale = 4;
                                double c_scale = System.Convert.ToDouble(780) / System.Convert.ToDouble(nc);
                                double r_scale = System.Convert.ToDouble(330) / System.Convert.ToDouble(nr);
                                if (c_scale < r_scale) { graphics_scale = System.Convert.ToInt32(Math.Floor(c_scale)); }
                                else { graphics_scale = System.Convert.ToInt32(Math.Floor(r_scale)); }
                                if (graphics_scale < 1) { graphics_scale = 1; }
                                m_objDrawingSurface = new Bitmap(nc * graphics_scale, nr * graphics_scale, System.Drawing.Imaging.PixelFormat.Format24bppRgb);
                                mygraphics = this.CreateGraphics();
                                Mapwindow.Visible = true;
                                //map_controls.Visible = true; */
                            }
                            catch { Debug.WriteLine("graphics initialisation failed "); input_data_error = true; }
                            if (input_data_error == false)
                            {
                                int count_intervene = 0;
                                t_intervene = 0;
                                if (t_intervene > 0) { read_soil_elevation_distance_from_output(t_intervene, workdir); }

                                //begining of looping
                                for (t = t_intervene; t < end_time; t++)
                                {

                                    try
                                    {
                                        every_timestep();
                                    }
                                    catch
                                    {
                                        Debug.WriteLine("failed to run in timestep " + t);
                                        // Catch for when the model crashes due to unknown reasons. The model will read the latest output and start calculating again from there which I named an intervention). When the crash occurs five times, the model breaks MM
                                        if (count_intervene < 5)
                                        {
                                            count_intervene += 1;
                                            t_intervene = t - (t % (int.Parse(Box_years_output.Text)));
                                            Debug.WriteLine("intervening at t" + t_intervene);
                                            read_soil_elevation_distance_from_output(t_intervene, workdir);
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
                    if (input_data_error == true)
                    {
                        MessageBox.Show("input data error - program can not yet run");
                        tabControl1.Visible = true;
                        Mapselector.Enabled = false;
                        map_controls.Visible = false;
                    }

                    if (Calibration_button.Checked == true)
                    {
                        //calculate how good this run was:
                        double current_error = calib_objective_function();
                        //store that information along with the parameter values used to achieve it:
                        calib_update_report(current_error);
                        if (current_error < best_error) { best_error = current_error; calib_update_best_paras(); best_run = run_number; }
                        //and check whether one 'level' of calibration has finished. If so, we have to change parameter values
                        Debug.WriteLine("run " + run_number + " number paras " + user_specified_number_of_calibration_parameters + " number ratios " + calibration_ratios_textbox.Text.Split(';').Length);
                        if ((run_number + 1) % Convert.ToInt32(Math.Pow(calibration_ratios_textbox.Text.Split(';').Length, user_specified_number_of_calibration_parameters)) == 0)
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
                                calib_shift_and_zoom(0, double.Parse(calibration_ratio_reduction_parameter_textbox.Text), double.Parse(parameter_K_textbox.Text));
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
                Debug.WriteLine("Error in accessing file " + this.dtm_input_filename_textbox.Text);
            }

        }  //end main

    }
}
