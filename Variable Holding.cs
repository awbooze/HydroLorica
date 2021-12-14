using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading;

namespace LORICAVariables
{
    public class GUIVariables
    {
        public GUIVariables(DelUpdateAllFields UAF, DelUpdateStatusPannel UPS, DelUpdateTimePannel UTP, DelUpdateVariables UF,
                            DelUpdateTimeSeries uts, DelUpdateProfile up, DelUpdateLanduse_determinator uld, DelUpdateSoildata usd,
                            DelDrawMap DM)
        {
            updateAllFields = UAF;
            updateStatusPannel = UPS;
            updateTimePannel = UTP;
            updateVariables = UF;

            Timeseries = new output_timeseries(uts);
            Profile = new output_profile(up);
            Landuse_determinator = new landuse_determinator(uld);
            Soildata = new soil_specifier(usd);

            drawMap = DM;

            //updateVariables();
        }

        public delegate void DelUpdateAllFields(); //Pushes Variables to GUI
        public delegate void DelUpdateStatusPannel(); //Updates GUI Status Pannel
        public delegate void DelUpdateTimePannel(); //Updates GUI Status Pannel
        public delegate void DelUpdateVariables(); //Pulls Values from GUI

        DelUpdateAllFields updateAllFields;
        DelUpdateStatusPannel updateStatusPannel;
        DelUpdateTimePannel updateTimePannel;
        DelUpdateVariables updateVariables;
        public void UpdateAllFields()
        {
            updateAllFields();
        }
        public void UpdateStatusPannel()
        {
            updateStatusPannel();
        }
        public void UpdateTimePannel()
        {
            updateTimePannel();
        }
        public void UpdateFields()
        {
            updateVariables();
        }


        public delegate void DelUpdateTimeSeries(); //Pushes Variables to TimeSeries
        public delegate void DelUpdateProfile(); //Pushes Variables to Profile
        public delegate void DelUpdateLanduse_determinator(); //Pushes Variables to Landuse_determinator
        public delegate void DelUpdateSoildata(); //Pushes Variables to Soildata


        //public delegate double Delbulk_Density_Calc(double, double , double, double, double, double, double, double); //Pulls Values from GUI
        //public Delbulk_Density_Calc bulk_density_calc;

        //if none of these fields get additional instructions, remove and rename the delegates
        //Ex: updateAllFields becomes UpdateAllFields




        public class output_timeseries
        {
            static DelUpdateTimeSeries updateTimeSeries;
            public static void UpdateTimeSeries()
            {
                updateTimeSeries();
            }
            public output_timeseries(DelUpdateTimeSeries uts)
            {
                updateTimeSeries = uts;
            }

            ReaderWriterLock Total_chem_weath_checkboxRWL = new ReaderWriterLock();
            protected bool total_chem_weath_checkbox = false;
            public bool Total_chem_weath_checkbox
            {
                get
                {
                    Total_chem_weath_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = total_chem_weath_checkbox;
                    Total_chem_weath_checkboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Total_chem_weath_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                    total_chem_weath_checkbox = value;
                    Total_chem_weath_checkboxRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }
            ReaderWriterLock Timeseries_number_soil_thicker_checkboxRWL = new ReaderWriterLock();
            protected bool timeseries_number_soil_thicker_checkbox = false;
            public bool Timeseries_number_soil_thicker_checkbox
            {
                get
                {
                    Timeseries_number_soil_thicker_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = timeseries_number_soil_thicker_checkbox;
                    Timeseries_number_soil_thicker_checkboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_number_soil_thicker_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_number_soil_thicker_checkbox = value;
                    Timeseries_number_soil_thicker_checkboxRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_soil_thicker_textboxRWL = new ReaderWriterLock();
            protected string timeseries_soil_thicker_textbox = "";
            public string Timeseries_soil_thicker_textbox
            {
                get
                {
                    Timeseries_soil_thicker_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = timeseries_soil_thicker_textbox;
                    Timeseries_soil_thicker_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_soil_thicker_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_soil_thicker_textbox = value;
                    Timeseries_soil_thicker_textboxRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Total_average_soilthickness_checkboxRWL = new ReaderWriterLock();
            protected bool total_average_soilthickness_checkbox = false;
            public bool Total_average_soilthickness_checkbox
            {
                get
                {
                    Total_average_soilthickness_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = total_average_soilthickness_checkbox;
                    Total_average_soilthickness_checkboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Total_average_soilthickness_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                    total_average_soilthickness_checkbox = value;
                    Total_average_soilthickness_checkboxRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_soil_depth_checkboxRWL = new ReaderWriterLock();
            protected bool timeseries_soil_depth_checkbox = false;
            public bool Timeseries_soil_depth_checkbox
            {
                get
                {
                    Timeseries_soil_depth_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = timeseries_soil_depth_checkbox;
                    Timeseries_soil_depth_checkboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_soil_depth_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_soil_depth_checkbox = value;
                    Timeseries_soil_depth_checkboxRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_soil_cell_rowRWL = new ReaderWriterLock();
            protected string timeseries_soil_cell_row = "";
            public string Timeseries_soil_cell_row
            {
                get
                {
                    Timeseries_soil_cell_rowRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = timeseries_soil_cell_row;
                    Timeseries_soil_cell_rowRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_soil_cell_rowRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_soil_cell_row = value;
                    Timeseries_soil_cell_rowRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_soil_cell_colRWL = new ReaderWriterLock();
            protected string timeseries_soil_cell_col = "";
            public string Timeseries_soil_cell_col
            {
                get
                {
                    Timeseries_soil_cell_colRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = timeseries_soil_cell_col;
                    Timeseries_soil_cell_colRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_soil_cell_colRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_soil_cell_col = value;
                    Timeseries_soil_cell_colRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_cell_waterflow_checkRWL = new ReaderWriterLock();
            protected bool timeseries_cell_waterflow_check = false;
            public bool Timeseries_cell_waterflow_check
            {
                get
                {
                    Timeseries_cell_waterflow_checkRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = timeseries_cell_waterflow_check;
                    Timeseries_cell_waterflow_checkRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_cell_waterflow_checkRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_cell_waterflow_check = value;
                    Timeseries_cell_waterflow_checkRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Total_fine_formed_checkboxRWL = new ReaderWriterLock();
            protected bool total_fine_formed_checkbox = false;
            public bool Total_fine_formed_checkbox
            {
                get
                {
                    Total_fine_formed_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = total_fine_formed_checkbox;
                    Total_fine_formed_checkboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Total_fine_formed_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                    total_fine_formed_checkbox = value;
                    Total_fine_formed_checkboxRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Total_mass_bioturbed_checkboxRWL = new ReaderWriterLock();
            protected bool total_mass_bioturbed_checkbox = false;
            public bool Total_mass_bioturbed_checkbox
            {
                get
                {
                    Total_mass_bioturbed_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = total_mass_bioturbed_checkbox;
                    Total_mass_bioturbed_checkboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Total_mass_bioturbed_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                    total_mass_bioturbed_checkbox = value;
                    Total_mass_bioturbed_checkboxRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Total_OM_input_checkboxRWL = new ReaderWriterLock();
            protected bool total_OM_input_checkbox = false;
            public bool Total_OM_input_checkbox
            {
                get
                {
                    Total_OM_input_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = total_OM_input_checkbox;
                    Total_OM_input_checkboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Total_OM_input_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                    total_OM_input_checkbox = value;
                    Total_OM_input_checkboxRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Total_fine_eluviated_checkboxRWL = new ReaderWriterLock();
            protected bool total_fine_eluviated_checkbox = false;
            public bool Total_fine_eluviated_checkbox
            {
                get
                {
                    Total_fine_eluviated_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = total_fine_eluviated_checkbox;
                    Total_fine_eluviated_checkboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Total_fine_eluviated_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                    total_fine_eluviated_checkbox = value;
                    Total_fine_eluviated_checkboxRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_erosion_thresholdRWL = new ReaderWriterLock();
            protected double timeseries_erosion_threshold = 0;
            public double Timeseries_erosion_threshold
            {
                get
                {
                    Timeseries_erosion_thresholdRWL.AcquireReaderLock(Timeout.Infinite);
                    double temp = timeseries_erosion_threshold;
                    Timeseries_erosion_thresholdRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_erosion_thresholdRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_erosion_threshold = value;
                    Timeseries_erosion_thresholdRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_deposition_thresholdRWL = new ReaderWriterLock();
            protected double timeseries_deposition_threshold = 0;
            public double Timeseries_deposition_threshold
            {
                get
                {
                    Timeseries_deposition_thresholdRWL.AcquireReaderLock(Timeout.Infinite);
                    double temp = timeseries_deposition_threshold;
                    Timeseries_deposition_thresholdRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_deposition_thresholdRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_deposition_threshold = value;
                    Timeseries_deposition_thresholdRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_waterflow_thresholdRWL = new ReaderWriterLock();
            protected double timeseries_waterflow_threshold = 0;
            public double Timeseries_waterflow_threshold
            {
                get
                {
                    Timeseries_waterflow_thresholdRWL.AcquireReaderLock(Timeout.Infinite);
                    double temp = timeseries_waterflow_threshold;
                    Timeseries_waterflow_thresholdRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_waterflow_thresholdRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_waterflow_threshold = value;
                    Timeseries_waterflow_thresholdRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_net_ero_checkRWL = new ReaderWriterLock();
            protected bool timeseries_net_ero_check = false;
            public bool Timeseries_net_ero_check
            {
                get
                {
                    Timeseries_net_ero_checkRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = timeseries_net_ero_check;
                    Timeseries_net_ero_checkRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_net_ero_checkRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_net_ero_check = value;
                    Timeseries_net_ero_checkRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_textbox_cell_rowRWL = new ReaderWriterLock();
            protected string timeseries_textbox_cell_row = "";
            public string Timeseries_textbox_cell_row
            {
                get
                {
                    Timeseries_textbox_cell_rowRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = timeseries_textbox_cell_row;
                    Timeseries_textbox_cell_rowRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_textbox_cell_rowRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_textbox_cell_row = value;
                    Timeseries_textbox_cell_rowRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_textbox_cell_colRWL = new ReaderWriterLock();
            protected string timeseries_textbox_cell_col = "";
            public string Timeseries_textbox_cell_col
            {
                get
                {
                    Timeseries_textbox_cell_colRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = timeseries_textbox_cell_col;
                    Timeseries_textbox_cell_colRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_textbox_cell_colRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_textbox_cell_col = value;
                    Timeseries_textbox_cell_colRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_number_dep_checkRWL = new ReaderWriterLock();
            protected bool timeseries_number_dep_check = false;
            public bool Timeseries_number_dep_check
            {
                get
                {
                    Timeseries_number_dep_checkRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = timeseries_number_dep_check;
                    Timeseries_number_dep_checkRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_number_dep_checkRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_number_dep_check = value;
                    Timeseries_number_dep_checkRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_cell_altitude_checkRWL = new ReaderWriterLock();
            protected bool timeseries_cell_altitude_check = false;
            public bool Timeseries_cell_altitude_check
            {
                get
                {
                    Timeseries_cell_altitude_checkRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = timeseries_cell_altitude_check;
                    Timeseries_cell_altitude_checkRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_cell_altitude_checkRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_cell_altitude_check = value;
                    Timeseries_cell_altitude_checkRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_number_erosion_checkRWL = new ReaderWriterLock();
            protected bool timeseries_number_erosion_check = false;
            public bool Timeseries_number_erosion_check
            {
                get
                {
                    Timeseries_number_erosion_checkRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = timeseries_number_erosion_check;
                    Timeseries_number_erosion_checkRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_number_erosion_checkRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_number_erosion_check = value;
                    Timeseries_number_erosion_checkRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_number_waterflow_checkRWL = new ReaderWriterLock();
            protected bool timeseries_number_waterflow_check = false;
            public bool Timeseries_number_waterflow_check
            {
                get
                {
                    Timeseries_number_waterflow_checkRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = timeseries_number_waterflow_check;
                    Timeseries_number_waterflow_checkRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_number_waterflow_checkRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_number_waterflow_check = value;
                    Timeseries_number_waterflow_checkRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_SDR_checkRWL = new ReaderWriterLock();
            protected bool timeseries_SDR_check = false;
            public bool Timeseries_SDR_check
            {
                get
                {
                    Timeseries_SDR_checkRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = timeseries_SDR_check;
                    Timeseries_SDR_checkRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_SDR_checkRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_SDR_check = value;
                    Timeseries_SDR_checkRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_total_average_alt_checkRWL = new ReaderWriterLock();
            protected bool timeseries_total_average_alt_check = false;
            public bool Timeseries_total_average_alt_check
            {
                get
                {
                    Timeseries_total_average_alt_checkRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = timeseries_total_average_alt_check;
                    Timeseries_total_average_alt_checkRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_total_average_alt_checkRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_total_average_alt_check = value;
                    Timeseries_total_average_alt_checkRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_total_dep_checkRWL = new ReaderWriterLock();
            protected bool timeseries_total_dep_check = false;
            public bool Timeseries_total_dep_check
            {
                get
                {
                    Timeseries_total_dep_checkRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = timeseries_total_dep_check;
                    Timeseries_total_dep_checkRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_total_dep_checkRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_total_dep_check = value;
                    Timeseries_total_dep_checkRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_total_ero_checkRWL = new ReaderWriterLock();
            protected bool timeseries_total_ero_check = false;
            public bool Timeseries_total_ero_check
            {
                get
                {
                    Timeseries_total_ero_checkRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = timeseries_total_ero_check;
                    Timeseries_total_ero_checkRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_total_ero_checkRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_total_ero_check = value;
                    Timeseries_total_ero_checkRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_total_evap_checkRWL = new ReaderWriterLock();
            protected bool timeseries_total_evap_check = false;
            public bool Timeseries_total_evap_check
            {
                get
                {
                    Timeseries_total_evap_checkRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = timeseries_total_evap_check;
                    Timeseries_total_evap_checkRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_total_evap_checkRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_total_evap_check = value;
                    Timeseries_total_evap_checkRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_total_infil_checkRWL = new ReaderWriterLock();
            protected bool timeseries_total_infil_check = false;
            public bool Timeseries_total_infil_check
            {
                get
                {
                    Timeseries_total_infil_checkRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = timeseries_total_infil_check;
                    Timeseries_total_infil_checkRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_total_infil_checkRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_total_infil_check = value;
                    Timeseries_total_infil_checkRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_total_outflow_checkRWL = new ReaderWriterLock();
            protected bool timeseries_total_outflow_check = false;
            public bool Timeseries_total_outflow_check
            {
                get
                {
                    Timeseries_total_outflow_checkRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = timeseries_total_outflow_check;
                    Timeseries_total_outflow_checkRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_total_outflow_checkRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_total_outflow_check = value;
                    Timeseries_total_outflow_checkRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_total_rain_checkRWL = new ReaderWriterLock();
            protected bool timeseries_total_rain_check = false;
            public bool Timeseries_total_rain_check
            {
                get
                {
                    Timeseries_total_rain_checkRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = timeseries_total_rain_check;
                    Timeseries_total_rain_checkRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_total_rain_checkRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_total_rain_check = value;
                    Timeseries_total_rain_checkRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Total_phys_weath_checkboxRWL = new ReaderWriterLock();
            protected bool total_phys_weath_checkbox = false;
            public bool Total_phys_weath_checkbox
            {
                get
                {
                    Total_phys_weath_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = total_phys_weath_checkbox;
                    Total_phys_weath_checkboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Total_phys_weath_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                    total_phys_weath_checkbox = value;
                    Total_phys_weath_checkboxRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_coarser_checkboxRWL = new ReaderWriterLock();
            protected bool timeseries_coarser_checkbox = false;
            public bool Timeseries_coarser_checkbox
            {
                get
                {
                    Timeseries_coarser_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = timeseries_coarser_checkbox;
                    Timeseries_coarser_checkboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_coarser_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_coarser_checkbox = value;
                    Timeseries_coarser_checkboxRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_soil_mass_checkboxRWL = new ReaderWriterLock();
            protected bool timeseries_soil_mass_checkbox = false;
            public bool Timeseries_soil_mass_checkbox
            {
                get
                {
                    Timeseries_soil_mass_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = timeseries_soil_mass_checkbox;
                    Timeseries_soil_mass_checkboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_soil_mass_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_soil_mass_checkbox = value;
                    Timeseries_soil_mass_checkboxRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

            ReaderWriterLock Timeseries_soil_coarser_fraction_textboxRWL = new ReaderWriterLock();
            protected string timeseries_soil_coarser_fraction_textbox = "";
            public string Timeseries_soil_coarser_fraction_textbox
            {
                get
                {
                    Timeseries_soil_coarser_fraction_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = timeseries_soil_coarser_fraction_textbox;
                    Timeseries_soil_coarser_fraction_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Timeseries_soil_coarser_fraction_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    timeseries_soil_coarser_fraction_textbox = value;
                    Timeseries_soil_coarser_fraction_textboxRWL.ReleaseWriterLock();

                    UpdateTimeSeries();
                }
            }

        }
        public class output_profile
        {
            static DelUpdateProfile updateProfile;
            public void UpdateProfile()
            {
                updateProfile();
            }
            public output_profile(DelUpdateProfile up)
            {
                updateProfile = up;
            }

            ReaderWriterLock Radio_pro1_colRWL = new ReaderWriterLock();
            protected bool radio_pro1_col = false;
            public bool Radio_pro1_col
            {
                get
                {
                    Radio_pro1_colRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = radio_pro1_col;
                    Radio_pro1_colRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Radio_pro1_colRWL.AcquireWriterLock(Timeout.Infinite);
                    radio_pro1_col = value;
                    Radio_pro1_colRWL.ReleaseWriterLock();

                    UpdateProfile();
                }
            }

            ReaderWriterLock Check_altitude_profile1RWL = new ReaderWriterLock();
            protected bool check_altitude_profile1 = false;
            public bool Check_altitude_profile1
            {
                get
                {
                    Check_altitude_profile1RWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = check_altitude_profile1;
                    Check_altitude_profile1RWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Check_altitude_profile1RWL.AcquireWriterLock(Timeout.Infinite);
                    check_altitude_profile1 = value;
                    Check_altitude_profile1RWL.ReleaseWriterLock();

                    UpdateProfile();
                }
            }

            ReaderWriterLock Check_waterflow_profile1RWL = new ReaderWriterLock();
            protected bool check_waterflow_profile1 = false;
            public bool Check_waterflow_profile1
            {
                get
                {
                    Check_waterflow_profile1RWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = check_waterflow_profile1;
                    Check_waterflow_profile1RWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Check_waterflow_profile1RWL.AcquireWriterLock(Timeout.Infinite);
                    check_waterflow_profile1 = value;
                    Check_waterflow_profile1RWL.ReleaseWriterLock();

                    UpdateProfile();
                }
            }

            ReaderWriterLock Radio_pro1_rowRWL = new ReaderWriterLock();
            protected bool radio_pro1_row = false;
            public bool Radio_pro1_row
            {
                get
                {
                    Radio_pro1_rowRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = radio_pro1_row;
                    Radio_pro1_rowRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Radio_pro1_rowRWL.AcquireWriterLock(Timeout.Infinite);
                    radio_pro1_row = value;
                    Radio_pro1_rowRWL.ReleaseWriterLock();

                    UpdateProfile();
                }
            }

            ReaderWriterLock Radio_pro2_colRWL = new ReaderWriterLock();
            protected bool radio_pro2_col = false;
            public bool Radio_pro2_col
            {
                get
                {
                    Radio_pro2_colRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = radio_pro2_col;
                    Radio_pro2_colRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Radio_pro2_colRWL.AcquireWriterLock(Timeout.Infinite);
                    radio_pro2_col = value;
                    Radio_pro2_colRWL.ReleaseWriterLock();

                    UpdateProfile();
                }
            }

            ReaderWriterLock Radio_pro2_rowRWL = new ReaderWriterLock();
            protected bool radio_pro2_row = false;
            public bool Radio_pro2_row
            {
                get
                {
                    Radio_pro2_rowRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = radio_pro2_row;
                    Radio_pro2_rowRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Radio_pro2_rowRWL.AcquireWriterLock(Timeout.Infinite);
                    radio_pro2_row = value;
                    Radio_pro2_rowRWL.ReleaseWriterLock();

                    UpdateProfile();
                }
            }

            ReaderWriterLock Radio_pro3_colRWL = new ReaderWriterLock();
            protected bool radio_pro3_col = false;
            public bool Radio_pro3_col
            {
                get
                {
                    Radio_pro3_colRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = radio_pro3_col;
                    Radio_pro3_colRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Radio_pro3_colRWL.AcquireWriterLock(Timeout.Infinite);
                    radio_pro3_col = value;
                    Radio_pro3_colRWL.ReleaseWriterLock();

                    UpdateProfile();
                }
            }

            ReaderWriterLock Radio_pro3_rowRWL = new ReaderWriterLock();
            protected bool radio_pro3_row = false;
            public bool Radio_pro3_row
            {
                get
                {
                    Radio_pro3_rowRWL.AcquireReaderLock(Timeout.Infinite);
                    bool temp = radio_pro3_row;
                    Radio_pro3_rowRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    Radio_pro3_rowRWL.AcquireWriterLock(Timeout.Infinite);
                    radio_pro3_row = value;
                    Radio_pro3_rowRWL.ReleaseWriterLock();

                    UpdateProfile();
                }
            }

            ReaderWriterLock P1_row_col_boxRWL = new ReaderWriterLock();
            protected string p1_row_col_box = "";
            public string P1_row_col_box
            {
                get
                {
                    P1_row_col_boxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = p1_row_col_box;
                    P1_row_col_boxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    P1_row_col_boxRWL.AcquireWriterLock(Timeout.Infinite);
                    p1_row_col_box = value;
                    P1_row_col_boxRWL.ReleaseWriterLock();

                    UpdateProfile();
                }
            }

            ReaderWriterLock P2_row_col_boxRWL = new ReaderWriterLock();
            protected string p2_row_col_box = "";
            public string P2_row_col_box
            {
                get
                {
                    P2_row_col_boxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = p2_row_col_box;
                    P2_row_col_boxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    P2_row_col_boxRWL.AcquireWriterLock(Timeout.Infinite);
                    p2_row_col_box = value;
                    P2_row_col_boxRWL.ReleaseWriterLock();

                    UpdateProfile();
                }
            }

            ReaderWriterLock P3_row_col_boxRWL = new ReaderWriterLock();
            protected string p3_row_col_box = "";
            public string P3_row_col_box
            {
                get
                {
                    P3_row_col_boxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = p3_row_col_box;
                    P3_row_col_boxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    P3_row_col_boxRWL.AcquireWriterLock(Timeout.Infinite);
                    p3_row_col_box = value;
                    P3_row_col_boxRWL.ReleaseWriterLock();

                    UpdateProfile();
                }
            }

        }
        public class landuse_determinator
        {
            static DelUpdateLanduse_determinator updateLanduse_determinator;
            public void UpdateLanduse_determinator()
            {
                updateLanduse_determinator();
            }
            public landuse_determinator(DelUpdateLanduse_determinator uld)
            {
                updateLanduse_determinator = uld;
            }

            #region LU1
            ReaderWriterLock LU1_Inf_textboxRWL = new ReaderWriterLock();
            protected string lU1_Inf_textbox;
            public string LU1_Inf_textbox
            {
                get
                {
                    LU1_Inf_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU1_Inf_textbox;
                    LU1_Inf_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU1_Inf_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU1_Inf_textbox = value;
                    LU1_Inf_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU1_Evap_textboxRWL = new ReaderWriterLock();
            protected string lU1_Evap_textbox;
            public string LU1_Evap_textbox
            {
                get
                {
                    LU1_Evap_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU1_Evap_textbox;
                    LU1_Evap_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU1_Evap_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU1_Evap_textbox = value;
                    LU1_Evap_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU1_Ero_textboxRWL = new ReaderWriterLock();
            protected string lU1_Ero_textbox;
            public string LU1_Ero_textbox
            {
                get
                {
                    LU1_Ero_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU1_Ero_textbox;
                    LU1_Ero_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU1_Ero_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU1_Ero_textbox = value;
                    LU1_Ero_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU1_Dep_textboxRWL = new ReaderWriterLock();
            protected string lU1_Dep_textbox;
            public string LU1_Dep_textbox
            {
                get
                {
                    LU1_Dep_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU1_Dep_textbox;
                    LU1_Dep_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU1_Dep_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU1_Dep_textbox = value;
                    LU1_Dep_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }

            #endregion

            #region LU2
            ReaderWriterLock LU2_Inf_textboxRWL = new ReaderWriterLock();
            protected string lU2_Inf_textbox;
            public string LU2_Inf_textbox
            {
                get
                {
                    LU2_Inf_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU2_Inf_textbox;
                    LU2_Inf_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU2_Inf_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU2_Inf_textbox = value;
                    LU2_Inf_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU2_Evap_textboxRWL = new ReaderWriterLock();
            protected string lU2_Evap_textbox;
            public string LU2_Evap_textbox
            {
                get
                {
                    LU2_Evap_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU2_Evap_textbox;
                    LU2_Evap_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU2_Evap_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU2_Evap_textbox = value;
                    LU2_Evap_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU2_Ero_textboxRWL = new ReaderWriterLock();
            protected string lU2_Ero_textbox;
            public string LU2_Ero_textbox
            {
                get
                {
                    LU2_Ero_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU2_Ero_textbox;
                    LU2_Ero_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU2_Ero_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU2_Ero_textbox = value;
                    LU2_Ero_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU2_Dep_textboxRWL = new ReaderWriterLock();
            protected string lU2_Dep_textbox;
            public string LU2_Dep_textbox
            {
                get
                {
                    LU2_Dep_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU2_Dep_textbox;
                    LU2_Dep_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU2_Dep_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU2_Dep_textbox = value;
                    LU2_Dep_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }

            #endregion

            #region LU3
            ReaderWriterLock LU3_Inf_textboxRWL = new ReaderWriterLock();
            protected string lU3_Inf_textbox;
            public string LU3_Inf_textbox
            {
                get
                {
                    LU3_Inf_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU3_Inf_textbox;
                    LU3_Inf_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU3_Inf_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU3_Inf_textbox = value;
                    LU3_Inf_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU3_Evap_textboxRWL = new ReaderWriterLock();
            protected string lU3_Evap_textbox;
            public string LU3_Evap_textbox
            {
                get
                {
                    LU3_Evap_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU3_Evap_textbox;
                    LU3_Evap_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU3_Evap_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU3_Evap_textbox = value;
                    LU3_Evap_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU3_Ero_textboxRWL = new ReaderWriterLock();
            protected string lU3_Ero_textbox;
            public string LU3_Ero_textbox
            {
                get
                {
                    LU3_Ero_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU3_Ero_textbox;
                    LU3_Ero_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU3_Ero_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU3_Ero_textbox = value;
                    LU3_Ero_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU3_Dep_textboxRWL = new ReaderWriterLock();
            protected string lU3_Dep_textbox;
            public string LU3_Dep_textbox
            {
                get
                {
                    LU3_Dep_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU3_Dep_textbox;
                    LU3_Dep_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU3_Dep_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU3_Dep_textbox = value;
                    LU3_Dep_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }

            #endregion

            #region LU4
            ReaderWriterLock LU4_Inf_textboxRWL = new ReaderWriterLock();
            protected string lU4_Inf_textbox;
            public string LU4_Inf_textbox
            {
                get
                {
                    LU4_Inf_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU4_Inf_textbox;
                    LU4_Inf_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU4_Inf_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU4_Inf_textbox = value;
                    LU4_Inf_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU4_Evap_textboxRWL = new ReaderWriterLock();
            protected string lU4_Evap_textbox;
            public string LU4_Evap_textbox
            {
                get
                {
                    LU4_Evap_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU4_Evap_textbox;
                    LU4_Evap_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU4_Evap_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU4_Evap_textbox = value;
                    LU4_Evap_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU4_Ero_textboxRWL = new ReaderWriterLock();
            protected string lU4_Ero_textbox;
            public string LU4_Ero_textbox
            {
                get
                {
                    LU4_Ero_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU4_Ero_textbox;
                    LU4_Ero_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU4_Ero_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU4_Ero_textbox = value;
                    LU4_Ero_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU4_Dep_textboxRWL = new ReaderWriterLock();
            protected string lU4_Dep_textbox;
            public string LU4_Dep_textbox
            {
                get
                {
                    LU4_Dep_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU4_Dep_textbox;
                    LU4_Dep_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU4_Dep_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU4_Dep_textbox = value;
                    LU4_Dep_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }

            #endregion

            #region LU5
            ReaderWriterLock LU5_Inf_textboxRWL = new ReaderWriterLock();
            protected string lU5_Inf_textbox;
            public string LU5_Inf_textbox
            {
                get
                {
                    LU5_Inf_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU5_Inf_textbox;
                    LU5_Inf_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU5_Inf_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU5_Inf_textbox = value;
                    LU5_Inf_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU5_Evap_textboxRWL = new ReaderWriterLock();
            protected string lU5_Evap_textbox;
            public string LU5_Evap_textbox
            {
                get
                {
                    LU5_Evap_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU5_Evap_textbox;
                    LU5_Evap_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU5_Evap_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU5_Evap_textbox = value;
                    LU5_Evap_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU5_Ero_textboxRWL = new ReaderWriterLock();
            protected string lU5_Ero_textbox;
            public string LU5_Ero_textbox
            {
                get
                {
                    LU5_Ero_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU5_Ero_textbox;
                    LU5_Ero_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU5_Ero_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU5_Ero_textbox = value;
                    LU5_Ero_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU5_Dep_textboxRWL = new ReaderWriterLock();
            protected string lU5_Dep_textbox;
            public string LU5_Dep_textbox
            {
                get
                {
                    LU5_Dep_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU5_Dep_textbox;
                    LU5_Dep_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU5_Dep_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU5_Dep_textbox = value;
                    LU5_Dep_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }

            #endregion

            #region LU6
            ReaderWriterLock LU6_Inf_textboxRWL = new ReaderWriterLock();
            protected string lU6_Inf_textbox;
            public string LU6_Inf_textbox
            {
                get
                {
                    LU6_Inf_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU6_Inf_textbox;
                    LU6_Inf_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU6_Inf_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU6_Inf_textbox = value;
                    LU6_Inf_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU6_Evap_textboxRWL = new ReaderWriterLock();
            protected string lU6_Evap_textbox;
            public string LU6_Evap_textbox
            {
                get
                {
                    LU6_Evap_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU6_Evap_textbox;
                    LU6_Evap_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU6_Evap_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU6_Evap_textbox = value;
                    LU6_Evap_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU6_Ero_textboxRWL = new ReaderWriterLock();
            protected string lU6_Ero_textbox;
            public string LU6_Ero_textbox
            {
                get
                {
                    LU6_Ero_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU6_Ero_textbox;
                    LU6_Ero_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU6_Ero_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU6_Ero_textbox = value;
                    LU6_Ero_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU6_Dep_textboxRWL = new ReaderWriterLock();
            protected string lU6_Dep_textbox;
            public string LU6_Dep_textbox
            {
                get
                {
                    LU6_Dep_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU6_Dep_textbox;
                    LU6_Dep_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU6_Dep_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU6_Dep_textbox = value;
                    LU6_Dep_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }

            #endregion

            #region LU7
            ReaderWriterLock LU7_Inf_textboxRWL = new ReaderWriterLock();
            protected string lU7_Inf_textbox;
            public string LU7_Inf_textbox
            {
                get
                {
                    LU7_Inf_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU7_Inf_textbox;
                    LU7_Inf_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU7_Inf_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU7_Inf_textbox = value;
                    LU7_Inf_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU7_Evap_textboxRWL = new ReaderWriterLock();
            protected string lU7_Evap_textbox;
            public string LU7_Evap_textbox
            {
                get
                {
                    LU7_Evap_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU7_Evap_textbox;
                    LU7_Evap_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU7_Evap_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU7_Evap_textbox = value;
                    LU7_Evap_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU7_Ero_textboxRWL = new ReaderWriterLock();
            protected string lU7_Ero_textbox;
            public string LU7_Ero_textbox
            {
                get
                {
                    LU7_Ero_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU7_Ero_textbox;
                    LU7_Ero_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU7_Ero_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU7_Ero_textbox = value;
                    LU7_Ero_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU7_Dep_textboxRWL = new ReaderWriterLock();
            protected string lU7_Dep_textbox;
            public string LU7_Dep_textbox
            {
                get
                {
                    LU7_Dep_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU7_Dep_textbox;
                    LU7_Dep_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU7_Dep_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU7_Dep_textbox = value;
                    LU7_Dep_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }

            #endregion

            #region LU8
            ReaderWriterLock LU8_Inf_textboxRWL = new ReaderWriterLock();
            protected string lU8_Inf_textbox;
            public string LU8_Inf_textbox
            {
                get
                {
                    LU8_Inf_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU8_Inf_textbox;
                    LU8_Inf_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU8_Inf_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU8_Inf_textbox = value;
                    LU8_Inf_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU8_Evap_textboxRWL = new ReaderWriterLock();
            protected string lU8_Evap_textbox;
            public string LU8_Evap_textbox
            {
                get
                {
                    LU8_Evap_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU8_Evap_textbox;
                    LU8_Evap_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU8_Evap_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU8_Evap_textbox = value;
                    LU8_Evap_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU8_Ero_textboxRWL = new ReaderWriterLock();
            protected string lU8_Ero_textbox;
            public string LU8_Ero_textbox
            {
                get
                {
                    LU8_Ero_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU8_Ero_textbox;
                    LU8_Ero_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU8_Ero_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU8_Ero_textbox = value;
                    LU8_Ero_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU8_Dep_textboxRWL = new ReaderWriterLock();
            protected string lU8_Dep_textbox;
            public string LU8_Dep_textbox
            {
                get
                {
                    LU8_Dep_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU8_Dep_textbox;
                    LU8_Dep_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU8_Dep_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU8_Dep_textbox = value;
                    LU8_Dep_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }

            #endregion

            #region LU9
            ReaderWriterLock LU9_Inf_textboxRWL = new ReaderWriterLock();
            protected string lU9_Inf_textbox;
            public string LU9_Inf_textbox
            {
                get
                {
                    LU9_Inf_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU9_Inf_textbox;
                    LU9_Inf_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU9_Inf_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU9_Inf_textbox = value;
                    LU9_Inf_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU9_Evap_textboxRWL = new ReaderWriterLock();
            protected string lU9_Evap_textbox;
            public string LU9_Evap_textbox
            {
                get
                {
                    LU9_Evap_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU9_Evap_textbox;
                    LU9_Evap_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU9_Evap_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU9_Evap_textbox = value;
                    LU9_Evap_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU9_Ero_textboxRWL = new ReaderWriterLock();
            protected string lU9_Ero_textbox;
            public string LU9_Ero_textbox
            {
                get
                {
                    LU9_Ero_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU9_Ero_textbox;
                    LU9_Ero_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU9_Ero_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU9_Ero_textbox = value;
                    LU9_Ero_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU9_Dep_textboxRWL = new ReaderWriterLock();
            protected string lU9_Dep_textbox;
            public string LU9_Dep_textbox
            {
                get
                {
                    LU9_Dep_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU9_Dep_textbox;
                    LU9_Dep_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU9_Dep_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU9_Dep_textbox = value;
                    LU9_Dep_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }

            #endregion

            #region LU10
            ReaderWriterLock LU10_Inf_textboxRWL = new ReaderWriterLock();
            protected string lU10_Inf_textbox;
            public string LU10_Inf_textbox
            {
                get
                {
                    LU10_Inf_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU10_Inf_textbox;
                    LU10_Inf_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU10_Inf_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU10_Inf_textbox = value;
                    LU10_Inf_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU10_Evap_textboxRWL = new ReaderWriterLock();
            protected string lU10_Evap_textbox;
            public string LU10_Evap_textbox
            {
                get
                {
                    LU10_Evap_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU10_Evap_textbox;
                    LU10_Evap_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU10_Evap_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU10_Evap_textbox = value;
                    LU10_Evap_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU10_Ero_textboxRWL = new ReaderWriterLock();
            protected string lU10_Ero_textbox;
            public string LU10_Ero_textbox
            {
                get
                {
                    LU10_Ero_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU10_Ero_textbox;
                    LU10_Ero_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU10_Ero_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU10_Ero_textbox = value;
                    LU10_Ero_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }
            ReaderWriterLock LU10_Dep_textboxRWL = new ReaderWriterLock();
            protected string lU10_Dep_textbox;
            public string LU10_Dep_textbox
            {
                get
                {
                    LU10_Dep_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = lU10_Dep_textbox;
                    LU10_Dep_textboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LU10_Dep_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                    lU10_Dep_textbox = value;
                    LU10_Dep_textboxRWL.ReleaseWriterLock();

                    UpdateLanduse_determinator();
                }
            }

            #endregion



        }
        public class soil_specifier
        {
            static DelUpdateSoildata updateSoildata;
            public void UpdateSoildata()
            {
                updateSoildata();
            }
            public soil_specifier(DelUpdateSoildata usd)
            {
                updateSoildata = usd;
            }

            ReaderWriterLock CoarseboxRWL = new ReaderWriterLock();
            protected string coarsebox;
            public string Coarsebox
            {
                get
                {
                    CoarseboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = coarsebox;
                    CoarseboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    CoarseboxRWL.AcquireWriterLock(Timeout.Infinite);
                    coarsebox = value;
                    CoarseboxRWL.ReleaseWriterLock();

                    UpdateSoildata();
                }
            }
            ReaderWriterLock SandboxRWL = new ReaderWriterLock();
            protected string sandbox;
            public string Sandbox
            {
                get
                {
                    SandboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = sandbox;
                    SandboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    SandboxRWL.AcquireWriterLock(Timeout.Infinite);
                    sandbox = value;
                    SandboxRWL.ReleaseWriterLock();

                    UpdateSoildata();
                }
            }
            ReaderWriterLock SiltboxRWL = new ReaderWriterLock();
            protected string siltbox;
            public string Siltbox
            {
                get
                {
                    SiltboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = siltbox;
                    SiltboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    SiltboxRWL.AcquireWriterLock(Timeout.Infinite);
                    siltbox = value;
                    SiltboxRWL.ReleaseWriterLock();

                    UpdateSoildata();
                }
            }
            ReaderWriterLock ClayboxRWL = new ReaderWriterLock();
            protected string claybox;
            public string Claybox
            {
                get
                {
                    ClayboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = claybox;
                    ClayboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    ClayboxRWL.AcquireWriterLock(Timeout.Infinite);
                    claybox = value;
                    ClayboxRWL.ReleaseWriterLock();

                    UpdateSoildata();
                }
            }
            ReaderWriterLock FineclayboxRWL = new ReaderWriterLock();
            protected string fineclaybox;
            public string Fineclaybox
            {
                get
                {
                    FineclayboxRWL.AcquireReaderLock(Timeout.Infinite);
                    string temp = fineclaybox;
                    FineclayboxRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    FineclayboxRWL.AcquireWriterLock(Timeout.Infinite);
                    fineclaybox = value;
                    FineclayboxRWL.ReleaseWriterLock();

                    UpdateSoildata();
                }
            }

        }




        ReaderWriterLock TimeseriesRWL = new ReaderWriterLock();
        protected output_timeseries timeseries;
        public output_timeseries Timeseries
        {
            get
            {
                TimeseriesRWL.AcquireReaderLock(Timeout.Infinite);
                output_timeseries temp = timeseries;
                TimeseriesRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                TimeseriesRWL.AcquireWriterLock(Timeout.Infinite);
                timeseries = value;
                TimeseriesRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock ProfileRWL = new ReaderWriterLock();
        protected output_profile profile;
        public output_profile Profile
        {
            get
            {
                ProfileRWL.AcquireReaderLock(Timeout.Infinite);
                output_profile temp = profile;
                ProfileRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                ProfileRWL.AcquireWriterLock(Timeout.Infinite);
                profile = value;
                ProfileRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Landuse_determinatorRWL = new ReaderWriterLock();
        protected landuse_determinator landuse;
        public landuse_determinator Landuse_determinator
        {
            get
            {
                Landuse_determinatorRWL.AcquireReaderLock(Timeout.Infinite);
                landuse_determinator temp = landuse;
                Landuse_determinatorRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Landuse_determinatorRWL.AcquireWriterLock(Timeout.Infinite);
                landuse = value;
                Landuse_determinatorRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock SoildataRWL = new ReaderWriterLock();
        protected soil_specifier soildata;
        public soil_specifier Soildata
        {
            get
            {
                SoildataRWL.AcquireReaderLock(Timeout.Infinite);
                soil_specifier temp = soildata;
                SoildataRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                SoildataRWL.AcquireWriterLock(Timeout.Infinite);
                soildata = value;
                SoildataRWL.ReleaseWriterLock();
            }
        }


        public class mapwindow
        {
            public class size
            {
                ReaderWriterLock WidthRWL = new ReaderWriterLock();
                protected int width;
                public int Width
                {
                    get
                    {
                        WidthRWL.AcquireReaderLock(Timeout.Infinite);
                        int temp = width;
                        WidthRWL.ReleaseReaderLock();

                        return temp;
                    }
                    set
                    {
                        WidthRWL.AcquireWriterLock(Timeout.Infinite);
                        width = value;
                        WidthRWL.ReleaseWriterLock();
                    }
                }
                ReaderWriterLock HeightRWL = new ReaderWriterLock();
                protected int height;
                public int Height
                {
                    get
                    {
                        HeightRWL.AcquireReaderLock(Timeout.Infinite);
                        int temp = height;
                        HeightRWL.ReleaseReaderLock();

                        return temp;
                    }
                    set
                    {
                        HeightRWL.AcquireWriterLock(Timeout.Infinite);
                        height = value;
                        HeightRWL.ReleaseWriterLock();
                    }
                }
            }
            public class location
            {
                ReaderWriterLock XRWL = new ReaderWriterLock();
                protected int x;
                public int X
                {
                    get
                    {
                        XRWL.AcquireReaderLock(Timeout.Infinite);
                        int temp = x;
                        XRWL.ReleaseReaderLock();

                        return temp;
                    }
                    set
                    {
                        XRWL.AcquireWriterLock(Timeout.Infinite);
                        x = value;
                        XRWL.ReleaseWriterLock();
                    }
                }
                ReaderWriterLock YRWL = new ReaderWriterLock();
                protected int y;
                public int Y
                {
                    get
                    {
                        YRWL.AcquireReaderLock(Timeout.Infinite);
                        int temp = y;
                        YRWL.ReleaseReaderLock();

                        return temp;
                    }
                    set
                    {
                        YRWL.AcquireWriterLock(Timeout.Infinite);
                        y = value;
                        YRWL.ReleaseWriterLock();
                    }
                }
            }

            ReaderWriterLock SizeRWL = new ReaderWriterLock();
            protected size s = new size();
            public size Size
            {
                get
                {
                    SizeRWL.AcquireReaderLock(Timeout.Infinite);
                    size temp = s;
                    SizeRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    SizeRWL.AcquireWriterLock(Timeout.Infinite);
                    s = value;
                    SizeRWL.ReleaseWriterLock();
                }
            }

            ReaderWriterLock LocationRWL = new ReaderWriterLock();
            protected location l = new location();
            public location Location
            {
                get
                {
                    LocationRWL.AcquireReaderLock(Timeout.Infinite);
                    location temp = l;
                    LocationRWL.ReleaseReaderLock();

                    return temp;
                }
                set
                {
                    LocationRWL.AcquireWriterLock(Timeout.Infinite);
                    l = value;
                    LocationRWL.ReleaseWriterLock();
                }
            }


        }

        ReaderWriterLock MapwindowRWL = new ReaderWriterLock();
        protected mapwindow mp = new mapwindow();
        public mapwindow Mapwindow
        {
            get
            {
                MapwindowRWL.AcquireReaderLock(Timeout.Infinite);
                mapwindow temp = mp;
                MapwindowRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                MapwindowRWL.AcquireWriterLock(Timeout.Infinite);
                mp = value;
                MapwindowRWL.ReleaseWriterLock();
            }
        }


        public delegate void DelDrawMap(Graphics graphics); //Pulls Values from GUI
        DelDrawMap drawMap;

        public void DrawMap(Graphics graphics)
        {
            UpdateFields();
            drawMap(graphics);
            updateVariables();
        }






        ReaderWriterLock End_timeRWL = new ReaderWriterLock();
        protected double end_time;
        public double End_time
        {
            get
            {
                End_timeRWL.AcquireReaderLock(Timeout.Infinite);
                double temp = end_time;
                End_timeRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                End_timeRWL.AcquireWriterLock(Timeout.Infinite);
                end_time = value;
                End_timeRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock DXRWL = new ReaderWriterLock();
        protected double dx;
        public double DX
        {
            get
            {
                DXRWL.AcquireReaderLock(Timeout.Infinite);
                double temp = dx;
                DXRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                DXRWL.AcquireWriterLock(Timeout.Infinite);
                dx = value;
                DXRWL.ReleaseWriterLock();
            }
        }






        ReaderWriterLock P_allRWL = new ReaderWriterLock();
        protected int[] p_all;
        public int[] P_all
        {
            get
            {
                P_allRWL.AcquireReaderLock(Timeout.Infinite);
                int[] temp = p_all;
                P_allRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                P_allRWL.AcquireWriterLock(Timeout.Infinite);
                p_all = value;
                P_allRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock ET0_allRWL = new ReaderWriterLock();
        protected int[] eT0_all;
        public int[] ET0_all
        {
            get
            {
                ET0_allRWL.AcquireReaderLock(Timeout.Infinite);
                int[] temp = eT0_all;
                ET0_allRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                ET0_allRWL.AcquireWriterLock(Timeout.Infinite);
                eT0_all = value;
                ET0_allRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock D_allRWL = new ReaderWriterLock();
        protected int[] d_all;
        public int[] D_all
        {
            get
            {
                D_allRWL.AcquireReaderLock(Timeout.Infinite);
                int[] temp = d_all;
                D_allRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                D_allRWL.AcquireWriterLock(Timeout.Infinite);
                d_all = value;
                D_allRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Tavg_allRWL = new ReaderWriterLock();
        protected int[] tavg_all;
        public int[] Tavg_all
        {
            get
            {
                Tavg_allRWL.AcquireReaderLock(Timeout.Infinite);
                int[] temp = tavg_all;
                Tavg_allRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Tavg_allRWL.AcquireWriterLock(Timeout.Infinite);
                tavg_all = value;
                Tavg_allRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Tmin_allRWL = new ReaderWriterLock();
        protected int[] tmin_all;
        public int[] Tmin_all
        {
            get
            {
                Tmin_allRWL.AcquireReaderLock(Timeout.Infinite);
                int[] temp = tmin_all;
                Tmin_allRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Tmin_allRWL.AcquireWriterLock(Timeout.Infinite);
                tmin_all = value;
                Tmin_allRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Tmax_allRWL = new ReaderWriterLock();
        protected int[] tmax_all;
        public int[] Tmax_all
        {
            get
            {
                Tmax_allRWL.AcquireReaderLock(Timeout.Infinite);
                int[] temp = tmax_all;
                Tmax_allRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Tmax_allRWL.AcquireWriterLock(Timeout.Infinite);
                tmax_all = value;
                Tmax_allRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock OFy_mRWL = new ReaderWriterLock();
        protected double[,,] oFy_m;
        public double[,,] OFy_m
        {
            get
            {
                OFy_mRWL.AcquireReaderLock(Timeout.Infinite);
                double[,,] temp = oFy_m;
                OFy_mRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                OFy_mRWL.AcquireWriterLock(Timeout.Infinite);
                oFy_m = value;
                OFy_mRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock CO3_kgRWL = new ReaderWriterLock();
        protected double[,,] cO3_kg;
        public double[,,] CO3_kg
        {
            get
            {
                CO3_kgRWL.AcquireReaderLock(Timeout.Infinite);
                double[,,] temp = cO3_kg;
                CO3_kgRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                CO3_kgRWL.AcquireWriterLock(Timeout.Infinite);
                cO3_kg = value;
                CO3_kgRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Rockweath_methodRWL = new ReaderWriterLock();
        protected int rockweath_method;
        public int Rockweath_method
        {
            get
            {
                Rockweath_methodRWL.AcquireReaderLock(Timeout.Infinite);
                int temp = rockweath_method;
                Rockweath_methodRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Rockweath_methodRWL.AcquireWriterLock(Timeout.Infinite);
                rockweath_method = value;
                Rockweath_methodRWL.ReleaseWriterLock();
            }
        }


        




















        //String Text Boxes


        ReaderWriterLock InfoStatusPanelRWL = new ReaderWriterLock();
        protected string infoStatusPanel = "";
        public string InfoStatusPanel
        {
            get
            {
                InfoStatusPanelRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = infoStatusPanel;
                InfoStatusPanelRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                InfoStatusPanelRWL.AcquireWriterLock(Timeout.Infinite);
                infoStatusPanel = value;
                InfoStatusPanelRWL.ReleaseWriterLock();

                UpdateStatusPannel();
            }
        }
        ReaderWriterLock TimeStatusPanelRWL = new ReaderWriterLock();
        protected string timeStatusPanel = "";
        public string TimeStatusPanel
        {
            get
            {
                TimeStatusPanelRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = timeStatusPanel;
                TimeStatusPanelRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                TimeStatusPanelRWL.AcquireWriterLock(Timeout.Infinite);
                timeStatusPanel = value;
                TimeStatusPanelRWL.ReleaseWriterLock();

                UpdateTimePannel();
            }
        }

        ReaderWriterLock Out_sed_statuspanelRWL = new ReaderWriterLock();
        protected string out_sed_statuspanel = "";
        public string Out_sed_statuspanel
        {
            get
            {
                Out_sed_statuspanelRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = out_sed_statuspanel;
                Out_sed_statuspanelRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Out_sed_statuspanelRWL.AcquireWriterLock(Timeout.Infinite);
                out_sed_statuspanel = value;
                Out_sed_statuspanelRWL.ReleaseWriterLock();

                UpdateStatusPannel();
            }
        }

        ReaderWriterLock DTM_input_filename_textboxRWL = new ReaderWriterLock();
        protected string dtm_input_filename_textbox = "";
        public string DTM_input_filename_textbox
        {
            get
            {
                DTM_input_filename_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = dtm_input_filename_textbox;
                DTM_input_filename_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                DTM_input_filename_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                dtm_input_filename_textbox = value;
                DTM_input_filename_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Number_runs_textboxRWL = new ReaderWriterLock();
        protected string number_runs_textbox = "";
        public string Number_runs_textbox
        {
            get
            {
                Number_runs_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = number_runs_textbox;
                Number_runs_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Number_runs_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                number_runs_textbox = value;
                Number_runs_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock ProcessStatusPanelRWL = new ReaderWriterLock();
        protected string processStatusPanel = "";
        public string ProcessStatusPanel
        {
            get
            {
                ProcessStatusPanelRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = processStatusPanel;
                ProcessStatusPanelRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                ProcessStatusPanelRWL.AcquireWriterLock(Timeout.Infinite);
                processStatusPanel = value;
                ProcessStatusPanelRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Calibration_ratios_textboxRWL = new ReaderWriterLock();
        protected string calibration_ratios_textbox = "";
        public string Calibration_ratios_textbox
        {
            get
            {
                Calibration_ratios_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = calibration_ratios_textbox;
                Calibration_ratios_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Calibration_ratios_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                calibration_ratios_textbox = value;
                Calibration_ratios_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock GoogAnimationSaveIntervalRWL = new ReaderWriterLock();
        protected string googAnimationSaveInterval = "";
        public string GoogAnimationSaveInterval
        {
            get
            {
                GoogAnimationSaveIntervalRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = googAnimationSaveInterval;
                GoogAnimationSaveIntervalRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                GoogAnimationSaveIntervalRWL.AcquireWriterLock(Timeout.Infinite);
                googAnimationSaveInterval = value;
                GoogAnimationSaveIntervalRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock UTMzoneboxRWL = new ReaderWriterLock();
        protected string uTMzonebox = "";
        public string UTMzonebox
        {
            get
            {
                UTMzoneboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = uTMzonebox;
                UTMzoneboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                UTMzoneboxRWL.AcquireWriterLock(Timeout.Infinite);
                uTMzonebox = value;
                UTMzoneboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock SaveintervalboxRWL = new ReaderWriterLock();
        protected string saveintervalbox = "";
        public string Saveintervalbox
        {
            get
            {
                SaveintervalboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = saveintervalbox;
                SaveintervalboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                SaveintervalboxRWL.AcquireWriterLock(Timeout.Infinite);
                saveintervalbox = value;
                SaveintervalboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Parameter_m_textboxRWL = new ReaderWriterLock();
        protected string parameter_m_textbox = "";
        public string Parameter_m_textbox
        {
            get
            {
                Parameter_m_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = parameter_m_textbox;
                Parameter_m_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Parameter_m_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                parameter_m_textbox = value;
                Parameter_m_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Parameter_n_textboxRWL = new ReaderWriterLock();
        protected string parameter_n_textbox = "";
        public string Parameter_n_textbox
        {
            get
            {
                Parameter_n_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = parameter_n_textbox;
                Parameter_n_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Parameter_n_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                parameter_n_textbox = value;
                Parameter_n_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Parameter_K_textboxRWL = new ReaderWriterLock();
        protected string parameter_K_textbox = "";
        public string Parameter_K_textbox
        {
            get
            {
                Parameter_K_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = parameter_K_textbox;
                Parameter_K_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Parameter_K_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                parameter_K_textbox = value;
                Parameter_K_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Bio_protection_constant_textboxRWL = new ReaderWriterLock();
        protected string bio_protection_constant_textbox = "";
        public string Bio_protection_constant_textbox
        {
            get
            {
                Bio_protection_constant_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = bio_protection_constant_textbox;
                Bio_protection_constant_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Bio_protection_constant_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                bio_protection_constant_textbox = value;
                Bio_protection_constant_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Rock_protection_constant_textboxRWL = new ReaderWriterLock();
        protected string rock_protection_constant_textbox = "";
        public string Rock_protection_constant_textbox
        {
            get
            {
                Rock_protection_constant_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = rock_protection_constant_textbox;
                Rock_protection_constant_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Rock_protection_constant_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                rock_protection_constant_textbox = value;
                Rock_protection_constant_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Selectivity_constant_textboxRWL = new ReaderWriterLock();
        protected string selectivity_constant_textbox = "";
        public string Selectivity_constant_textbox
        {
            get
            {
                Selectivity_constant_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = selectivity_constant_textbox;
                Selectivity_constant_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Selectivity_constant_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                selectivity_constant_textbox = value;
                Selectivity_constant_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Erosion_threshold_textboxRWL = new ReaderWriterLock();
        protected string erosion_threshold_textbox = "";
        public string Erosion_threshold_textbox
        {
            get
            {
                Erosion_threshold_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = erosion_threshold_textbox;
                Erosion_threshold_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Erosion_threshold_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                erosion_threshold_textbox = value;
                Erosion_threshold_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Parameter_ploughing_depth_textboxRWL = new ReaderWriterLock();
        protected string parameter_ploughing_depth_textbox = "";
        public string Parameter_ploughing_depth_textbox
        {
            get
            {
                Parameter_ploughing_depth_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = parameter_ploughing_depth_textbox;
                Parameter_ploughing_depth_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Parameter_ploughing_depth_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                parameter_ploughing_depth_textbox = value;
                Parameter_ploughing_depth_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Parameter_tillage_constant_textboxRWL = new ReaderWriterLock();
        protected string parameter_tillage_constant_textbox = "";
        public string Parameter_tillage_constant_textbox
        {
            get
            {
                Parameter_tillage_constant_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = parameter_tillage_constant_textbox;
                Parameter_tillage_constant_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Parameter_tillage_constant_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                parameter_tillage_constant_textbox = value;
                Parameter_tillage_constant_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Parameter_diffusivity_textboxRWL = new ReaderWriterLock();
        protected string parameter_diffusivity_textbox = "";
        public string Parameter_diffusivity_textbox
        {
            get
            {
                Parameter_diffusivity_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = parameter_diffusivity_textbox;
                Parameter_diffusivity_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Parameter_diffusivity_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                parameter_diffusivity_textbox = value;
                Parameter_diffusivity_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Parameter_P0_textboxRWL = new ReaderWriterLock();
        protected string parameter_P0_textbox = "";
        public string Parameter_P0_textbox
        {
            get
            {
                Parameter_P0_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = parameter_P0_textbox;
                Parameter_P0_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Parameter_P0_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                parameter_P0_textbox = value;
                Parameter_P0_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Parameter_k1_textboxRWL = new ReaderWriterLock();
        protected string parameter_k1_textbox = "";
        public string Parameter_k1_textbox
        {
            get
            {
                Parameter_k1_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = parameter_k1_textbox;
                Parameter_k1_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Parameter_k1_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                parameter_k1_textbox = value;
                Parameter_k1_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Parameter_k2_textboxRWL = new ReaderWriterLock();
        protected string parameter_k2_textbox = "";
        public string Parameter_k2_textbox
        {
            get
            {
                Parameter_k2_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = parameter_k2_textbox;
                Parameter_k2_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Parameter_k2_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                parameter_k2_textbox = value;
                Parameter_k2_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Parameter_Pa_textboxRWL = new ReaderWriterLock();
        protected string parameter_Pa_textbox = "";
        public string Parameter_Pa_textbox
        {
            get
            {
                Parameter_Pa_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = parameter_Pa_textbox;
                Parameter_Pa_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Parameter_Pa_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                parameter_Pa_textbox = value;
                Parameter_Pa_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Tilting_rate_textboxRWL = new ReaderWriterLock();
        protected string tilting_rate_textbox = "";
        public string Tilting_rate_textbox
        {
            get
            {
                Tilting_rate_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = tilting_rate_textbox;
                Tilting_rate_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Tilting_rate_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                tilting_rate_textbox = value;
                Tilting_rate_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Text_lift_row_lessRWL = new ReaderWriterLock();
        protected string text_lift_row_less = "";
        public string Text_lift_row_less
        {
            get
            {
                Text_lift_row_lessRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = text_lift_row_less;
                Text_lift_row_lessRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Text_lift_row_lessRWL.AcquireWriterLock(Timeout.Infinite);
                text_lift_row_less = value;
                Text_lift_row_lessRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Text_lift_row_moreRWL = new ReaderWriterLock();
        protected string text_lift_row_more = "";
        public string Text_lift_row_more
        {
            get
            {
                Text_lift_row_moreRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = text_lift_row_more;
                Text_lift_row_moreRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Text_lift_row_moreRWL.AcquireWriterLock(Timeout.Infinite);
                text_lift_row_more = value;
                Text_lift_row_moreRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Text_lift_col_lessRWL = new ReaderWriterLock();
        protected string text_lift_col_less = "";
        public string Text_lift_col_less
        {
            get
            {
                Text_lift_col_lessRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = text_lift_col_less;
                Text_lift_col_lessRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Text_lift_col_lessRWL.AcquireWriterLock(Timeout.Infinite);
                text_lift_col_less = value;
                Text_lift_col_lessRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Text_lift_col_moreRWL = new ReaderWriterLock();
        protected string text_lift_col_more = "";
        public string Text_lift_col_more
        {
            get
            {
                Text_lift_col_moreRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = text_lift_col_more;
                Text_lift_col_moreRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Text_lift_col_moreRWL.AcquireWriterLock(Timeout.Infinite);
                text_lift_col_more = value;
                Text_lift_col_moreRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Uplift_rate_textboxRWL = new ReaderWriterLock();
        protected string uplift_rate_textbox = "";
        public string Uplift_rate_textbox
        {
            get
            {
                Uplift_rate_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = uplift_rate_textbox;
                Uplift_rate_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Uplift_rate_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                uplift_rate_textbox = value;
                Uplift_rate_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Tf_WRWL = new ReaderWriterLock();
        protected string tf_W = "";
        public string Tf_W
        {
            get
            {
                Tf_WRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = tf_W;
                Tf_WRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Tf_WRWL.AcquireWriterLock(Timeout.Infinite);
                tf_W = value;
                Tf_WRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Tf_DRWL = new ReaderWriterLock();
        protected string tf_D = "";
        public string Tf_D
        {
            get
            {
                Tf_DRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = tf_D;
                Tf_DRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Tf_DRWL.AcquireWriterLock(Timeout.Infinite);
                tf_D = value;
                Tf_DRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Tf_growthRWL = new ReaderWriterLock();
        protected string tf_growth = "";
        public string Tf_growth
        {
            get
            {
                Tf_growthRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = tf_growth;
                Tf_growthRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Tf_growthRWL.AcquireWriterLock(Timeout.Infinite);
                tf_growth = value;
                Tf_growthRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Tf_ageRWL = new ReaderWriterLock();
        protected string tf_age = "";
        public string Tf_age
        {
            get
            {
                Tf_ageRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = tf_age;
                Tf_ageRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Tf_ageRWL.AcquireWriterLock(Timeout.Infinite);
                tf_age = value;
                Tf_ageRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Tf_freqRWL = new ReaderWriterLock();
        protected string tf_freq = "";
        public string Tf_freq
        {
            get
            {
                Tf_freqRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = tf_freq;
                Tf_freqRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Tf_freqRWL.AcquireWriterLock(Timeout.Infinite);
                tf_freq = value;
                Tf_freqRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Physical_weath_C1_textboxRWL = new ReaderWriterLock();
        protected string physical_weath_C1_textbox = "";
        public string Physical_weath_C1_textbox
        {
            get
            {
                Physical_weath_C1_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = physical_weath_C1_textbox;
                Physical_weath_C1_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Physical_weath_C1_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                physical_weath_C1_textbox = value;
                Physical_weath_C1_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Physical_weath_constant1RWL = new ReaderWriterLock();
        protected string physical_weath_constant1 = "";
        public string Physical_weath_constant1
        {
            get
            {
                Physical_weath_constant1RWL.AcquireReaderLock(Timeout.Infinite);
                string temp = physical_weath_constant1;
                Physical_weath_constant1RWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Physical_weath_constant1RWL.AcquireWriterLock(Timeout.Infinite);
                physical_weath_constant1 = value;
                Physical_weath_constant1RWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Physical_weath_constant2RWL = new ReaderWriterLock();
        protected string physical_weath_constant2 = "";
        public string Physical_weath_constant2
        {
            get
            {
                Physical_weath_constant2RWL.AcquireReaderLock(Timeout.Infinite);
                string temp = physical_weath_constant2;
                Physical_weath_constant2RWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Physical_weath_constant2RWL.AcquireWriterLock(Timeout.Infinite);
                physical_weath_constant2 = value;
                Physical_weath_constant2RWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Chem_weath_rate_constant_textboxRWL = new ReaderWriterLock();
        protected string chem_weath_rate_constant_textbox = "";
        public string Chem_weath_rate_constant_textbox
        {
            get
            {
                Chem_weath_rate_constant_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = chem_weath_rate_constant_textbox;
                Chem_weath_rate_constant_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Chem_weath_rate_constant_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                chem_weath_rate_constant_textbox = value;
                Chem_weath_rate_constant_textboxRWL.ReleaseWriterLock();
            }
        }


        ReaderWriterLock Chem_weath_depth_constant_textboxRWL = new ReaderWriterLock();
        protected string chem_weath_depth_constant_textbox = "";
        public string Chem_weath_depth_constant_textbox
        {
            get
            {
                Chem_weath_depth_constant_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = chem_weath_depth_constant_textbox;
                Chem_weath_depth_constant_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Chem_weath_depth_constant_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                chem_weath_depth_constant_textbox = value;
                Chem_weath_depth_constant_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Chem_weath_specific_coefficient_textboxRWL = new ReaderWriterLock();
        protected string chem_weath_specific_coefficient_textbox = "";
        public string Chem_weath_specific_coefficient_textbox
        {
            get
            {
                Chem_weath_specific_coefficient_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = chem_weath_specific_coefficient_textbox;
                Chem_weath_specific_coefficient_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Chem_weath_specific_coefficient_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                chem_weath_specific_coefficient_textbox = value;
                Chem_weath_specific_coefficient_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Specific_area_coarse_textboxRWL = new ReaderWriterLock();
        protected string specific_area_coarse_textbox = "";
        public string Specific_area_coarse_textbox
        {
            get
            {
                Specific_area_coarse_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = specific_area_coarse_textbox;
                Specific_area_coarse_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Specific_area_coarse_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                specific_area_coarse_textbox = value;
                Specific_area_coarse_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Specific_area_sand_textboxRWL = new ReaderWriterLock();
        protected string specific_area_sand_textbox = "";
        public string Specific_area_sand_textbox
        {
            get
            {
                Specific_area_sand_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = specific_area_sand_textbox;
                Specific_area_sand_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Specific_area_sand_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                specific_area_sand_textbox = value;
                Specific_area_sand_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Specific_area_silt_textboxRWL = new ReaderWriterLock();
        protected string specific_area_silt_textbox = "";
        public string Specific_area_silt_textbox
        {
            get
            {
                Specific_area_silt_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = specific_area_silt_textbox;
                Specific_area_silt_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Specific_area_silt_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                specific_area_silt_textbox = value;
                Specific_area_silt_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Specific_area_clay_textboxRWL = new ReaderWriterLock();
        protected string specific_area_clay_textbox = "";
        public string Specific_area_clay_textbox
        {
            get
            {
                Specific_area_clay_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = specific_area_clay_textbox;
                Specific_area_clay_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Specific_area_clay_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                specific_area_clay_textbox = value;
                Specific_area_clay_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Specific_area_fine_clay_textboxRWL = new ReaderWriterLock();
        protected string specific_area_fine_clay_textbox = "";
        public string Specific_area_fine_clay_textbox
        {
            get
            {
                Specific_area_fine_clay_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = specific_area_fine_clay_textbox;
                Specific_area_fine_clay_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Specific_area_fine_clay_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                specific_area_fine_clay_textbox = value;
                Specific_area_fine_clay_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Clay_neoform_constant_textboxRWL = new ReaderWriterLock();
        protected string clay_neoform_constant_textbox = "";
        public string Clay_neoform_constant_textbox
        {
            get
            {
                Clay_neoform_constant_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = clay_neoform_constant_textbox;
                Clay_neoform_constant_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Clay_neoform_constant_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                clay_neoform_constant_textbox = value;
                Clay_neoform_constant_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Clay_neoform_C1_textboxRWL = new ReaderWriterLock();
        protected string clay_neoform_C1_textbox = "";
        public string Clay_neoform_C1_textbox
        {
            get
            {
                Clay_neoform_C1_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = clay_neoform_C1_textbox;
                Clay_neoform_C1_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Clay_neoform_C1_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                clay_neoform_C1_textbox = value;
                Clay_neoform_C1_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Clay_neoform_C2_textboxRWL = new ReaderWriterLock();
        protected string clay_neoform_C2_textbox = "";
        public string Clay_neoform_C2_textbox
        {
            get
            {
                Clay_neoform_C2_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = clay_neoform_C2_textbox;
                Clay_neoform_C2_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Clay_neoform_C2_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                clay_neoform_C2_textbox = value;
                Clay_neoform_C2_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Maximum_eluviation_textboxRWL = new ReaderWriterLock();
        protected string maximum_eluviation_textbox = "";
        public string Maximum_eluviation_textbox
        {
            get
            {
                Maximum_eluviation_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = maximum_eluviation_textbox;
                Maximum_eluviation_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Maximum_eluviation_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                maximum_eluviation_textbox = value;
                Maximum_eluviation_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Eluviation_coefficient_textboxRWL = new ReaderWriterLock();
        protected string eluviation_coefficient_textbox = "";
        public string Eluviation_coefficient_textbox
        {
            get
            {
                Eluviation_coefficient_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = eluviation_coefficient_textbox;
                Eluviation_coefficient_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Eluviation_coefficient_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                eluviation_coefficient_textbox = value;
                Eluviation_coefficient_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Ct_depth_decayRWL = new ReaderWriterLock();
        protected string ct_depth_decay = "";
        public string Ct_depth_decay
        {
            get
            {
                Ct_depth_decayRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = ct_depth_decay;
                Ct_depth_decayRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Ct_depth_decayRWL.AcquireWriterLock(Timeout.Infinite);
                ct_depth_decay = value;
                Ct_depth_decayRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Potential_bioturbation_textboxRWL = new ReaderWriterLock();
        protected string potential_bioturbation_textbox = "";
        public string Potential_bioturbation_textbox
        {
            get
            {
                Potential_bioturbation_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = potential_bioturbation_textbox;
                Potential_bioturbation_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Potential_bioturbation_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                potential_bioturbation_textbox = value;
                Potential_bioturbation_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Bioturbation_depth_decay_textboxRWL = new ReaderWriterLock();
        protected string bioturbation_depth_decay_textbox = "";
        public string Bioturbation_depth_decay_textbox
        {
            get
            {
                Bioturbation_depth_decay_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = bioturbation_depth_decay_textbox;
                Bioturbation_depth_decay_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Bioturbation_depth_decay_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                bioturbation_depth_decay_textbox = value;
                Bioturbation_depth_decay_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Carbon_input_textboxRWL = new ReaderWriterLock();
        protected string carbon_input_textbox = "";
        public string Carbon_input_textbox
        {
            get
            {
                Carbon_input_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = carbon_input_textbox;
                Carbon_input_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Carbon_input_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                carbon_input_textbox = value;
                Carbon_input_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Carbon_depth_decay_textboxRWL = new ReaderWriterLock();
        protected string carbon_depth_decay_textbox = "";
        public string Carbon_depth_decay_textbox
        {
            get
            {
                Carbon_depth_decay_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = carbon_depth_decay_textbox;
                Carbon_depth_decay_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Carbon_depth_decay_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                carbon_depth_decay_textbox = value;
                Carbon_depth_decay_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Carbon_humification_fraction_textboxRWL = new ReaderWriterLock();
        protected string carbon_humification_fraction_textbox = "";
        public string Carbon_humification_fraction_textbox
        {
            get
            {
                Carbon_humification_fraction_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = carbon_humification_fraction_textbox;
                Carbon_humification_fraction_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Carbon_humification_fraction_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                carbon_humification_fraction_textbox = value;
                Carbon_humification_fraction_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Carbon_y_decomp_rate_textboxRWL = new ReaderWriterLock();
        protected string carbon_y_decomp_rate_textbox = "";
        public string Carbon_y_decomp_rate_textbox
        {
            get
            {
                Carbon_y_decomp_rate_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = carbon_y_decomp_rate_textbox;
                Carbon_y_decomp_rate_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Carbon_y_decomp_rate_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                carbon_y_decomp_rate_textbox = value;
                Carbon_y_decomp_rate_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Carbon_o_decomp_rate_textboxRWL = new ReaderWriterLock();
        protected string carbon_o_decomp_rate_textbox = "";
        public string Carbon_o_decomp_rate_textbox
        {
            get
            {
                Carbon_o_decomp_rate_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = carbon_o_decomp_rate_textbox;
                Carbon_o_decomp_rate_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Carbon_o_decomp_rate_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                carbon_o_decomp_rate_textbox = value;
                Carbon_o_decomp_rate_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Carbon_y_depth_decay_textboxRWL = new ReaderWriterLock();
        protected string carbon_y_depth_decay_textbox = "";
        public string Carbon_y_depth_decay_textbox
        {
            get
            {
                Carbon_y_depth_decay_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = carbon_y_depth_decay_textbox;
                Carbon_y_depth_decay_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Carbon_y_depth_decay_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                carbon_y_depth_decay_textbox = value;
                Carbon_y_depth_decay_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Carbon_o_twi_decay_textboxRWL = new ReaderWriterLock();
        protected string carbon_o_twi_decay_textbox = "";
        public string Carbon_o_twi_decay_textbox
        {
            get
            {
                Carbon_o_twi_decay_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = carbon_o_twi_decay_textbox;
                Carbon_o_twi_decay_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Carbon_o_twi_decay_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                carbon_o_twi_decay_textbox = value;
                Carbon_o_twi_decay_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Carbon_o_depth_decay_textboxRWL = new ReaderWriterLock();
        protected string carbon_o_depth_decay_textbox = "";
        public string Carbon_o_depth_decay_textbox
        {
            get
            {
                Carbon_o_depth_decay_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = carbon_o_depth_decay_textbox;
                Carbon_o_depth_decay_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Carbon_o_depth_decay_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                carbon_o_depth_decay_textbox = value;
                Carbon_o_depth_decay_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Carbon_y_twi_decay_textboxRWL = new ReaderWriterLock();
        protected string carbon_y_twi_decay_textbox = "";
        public string Carbon_y_twi_decay_textbox
        {
            get
            {
                Carbon_y_twi_decay_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = carbon_y_twi_decay_textbox;
                Carbon_y_twi_decay_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Carbon_y_twi_decay_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                carbon_y_twi_decay_textbox = value;
                Carbon_y_twi_decay_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Landuse_constant_value_boxRWL = new ReaderWriterLock();
        protected string landuse_constant_value_box = "";
        public string Landuse_constant_value_box
        {
            get
            {
                Landuse_constant_value_boxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = landuse_constant_value_box;
                Landuse_constant_value_boxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Landuse_constant_value_boxRWL.AcquireWriterLock(Timeout.Infinite);
                landuse_constant_value_box = value;
                Landuse_constant_value_boxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Evap_constant_value_boxRWL = new ReaderWriterLock();
        protected string evap_constant_value_box = "";
        public string Evap_constant_value_box
        {
            get
            {
                Evap_constant_value_boxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = evap_constant_value_box;
                Evap_constant_value_boxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Evap_constant_value_boxRWL.AcquireWriterLock(Timeout.Infinite);
                evap_constant_value_box = value;
                Evap_constant_value_boxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Infil_constant_value_boxRWL = new ReaderWriterLock();
        protected string infil_constant_value_box = "";
        public string Infil_constant_value_box
        {
            get
            {
                Infil_constant_value_boxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = infil_constant_value_box;
                Infil_constant_value_boxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Infil_constant_value_boxRWL.AcquireWriterLock(Timeout.Infinite);
                infil_constant_value_box = value;
                Infil_constant_value_boxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Rainfall_constant_value_boxRWL = new ReaderWriterLock();
        protected string rainfall_constant_value_box = "";
        public string Rainfall_constant_value_box
        {
            get
            {
                Rainfall_constant_value_boxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = rainfall_constant_value_box;
                Rainfall_constant_value_boxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Rainfall_constant_value_boxRWL.AcquireWriterLock(Timeout.Infinite);
                rainfall_constant_value_box = value;
                Rainfall_constant_value_boxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Temp_constant_value_boxRWL = new ReaderWriterLock();
        protected string temp_constant_value_box = "";
        public string Temp_constant_value_box
        {
            get
            {
                Temp_constant_value_boxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = temp_constant_value_box;
                Temp_constant_value_boxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Temp_constant_value_boxRWL.AcquireWriterLock(Timeout.Infinite);
                temp_constant_value_box = value;
                Temp_constant_value_boxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Calibration_ratio_reduction_parameter_textboxRWL = new ReaderWriterLock();
        protected string calibration_ratio_reduction_parameter_textbox = "";
        public string Calibration_ratio_reduction_parameter_textbox
        {
            get
            {
                Calibration_ratio_reduction_parameter_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = calibration_ratio_reduction_parameter_textbox;
                Calibration_ratio_reduction_parameter_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Calibration_ratio_reduction_parameter_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                calibration_ratio_reduction_parameter_textbox = value;
                Calibration_ratio_reduction_parameter_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Soildepth_constant_value_boxRWL = new ReaderWriterLock();
        protected string soildepth_constant_value_box = "";
        public string Soildepth_constant_value_box
        {
            get
            {
                Soildepth_constant_value_boxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = soildepth_constant_value_box;
                Soildepth_constant_value_boxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Soildepth_constant_value_boxRWL.AcquireWriterLock(Timeout.Infinite);
                soildepth_constant_value_box = value;
                Soildepth_constant_value_boxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Box_years_outputRWL = new ReaderWriterLock();
        protected string box_years_output = "";
        public string Box_years_output
        {
            get
            {
                Box_years_outputRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = box_years_output;
                Box_years_outputRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Box_years_outputRWL.AcquireWriterLock(Timeout.Infinite);
                box_years_output = value;
                Box_years_outputRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Ct_v0_JagercikovaRWL = new ReaderWriterLock();
        protected string ct_v0_Jagercikova = "";
        public string Ct_v0_Jagercikova
        {
            get
            {
                Ct_v0_JagercikovaRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = ct_v0_Jagercikova;
                Ct_v0_JagercikovaRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Ct_v0_JagercikovaRWL.AcquireWriterLock(Timeout.Infinite);
                ct_v0_Jagercikova = value;
                Ct_v0_JagercikovaRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Ct_dd_JagercikovaRWL = new ReaderWriterLock();
        protected string ct_dd_Jagercikova = "";
        public string Ct_dd_Jagercikova
        {
            get
            {
                Ct_dd_JagercikovaRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = ct_dd_Jagercikova;
                Ct_dd_JagercikovaRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Ct_dd_JagercikovaRWL.AcquireWriterLock(Timeout.Infinite);
                ct_dd_Jagercikova = value;
                Ct_dd_JagercikovaRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Upper_particle_coarse_textboxRWL = new ReaderWriterLock();
        protected string upper_particle_coarse_textbox = "";
        public string Upper_particle_coarse_textbox
        {
            get
            {
                Upper_particle_coarse_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = upper_particle_coarse_textbox;
                Upper_particle_coarse_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Upper_particle_coarse_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                Upper_particle_coarse_textbox = value;
                Upper_particle_coarse_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Upper_particle_sand_textboxRWL = new ReaderWriterLock();
        protected string upper_particle_sand_textbox = "";
        public string Upper_particle_sand_textbox
        {
            get
            {
                Upper_particle_sand_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = upper_particle_sand_textbox;
                Upper_particle_sand_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Upper_particle_sand_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                upper_particle_sand_textbox = value;
                Upper_particle_sand_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Upper_particle_silt_textboxRWL = new ReaderWriterLock();
        protected string upper_particle_silt_textbox = "";
        public string Upper_particle_silt_textbox
        {
            get
            {
                Upper_particle_silt_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = upper_particle_silt_textbox;
                Upper_particle_silt_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Upper_particle_silt_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                upper_particle_silt_textbox = value;
                Upper_particle_silt_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Upper_particle_clay_textboxRWL = new ReaderWriterLock();
        protected string upper_particle_clay_textbox = "";
        public string Upper_particle_clay_textbox
        {
            get
            {
                Upper_particle_clay_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = upper_particle_clay_textbox;
                Upper_particle_clay_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Upper_particle_clay_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                upper_particle_clay_textbox = value;
                Upper_particle_clay_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Upper_particle_fine_clay_textboxRWL = new ReaderWriterLock();
        protected string upper_particle_fine_clay_textbox = "";
        public string Upper_particle_fine_clay_textbox
        {
            get
            {
                Upper_particle_fine_clay_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = upper_particle_fine_clay_textbox;
                Upper_particle_fine_clay_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Upper_particle_fine_clay_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                upper_particle_fine_clay_textbox = value;
                Upper_particle_fine_clay_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Soildepth_input_filename_textboxRWL = new ReaderWriterLock();
        protected string soildepth_input_filename_textbox = "";
        public string Soildepth_input_filename_textbox
        {
            get
            {
                Soildepth_input_filename_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = soildepth_input_filename_textbox;
                Soildepth_input_filename_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Soildepth_input_filename_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                soildepth_input_filename_textbox = value;
                Soildepth_input_filename_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Snowmelt_factor_textboxRWL = new ReaderWriterLock();
        protected string snowmelt_factor_textbox = "";
        public string Snowmelt_factor_textbox
        {
            get
            {
                Snowmelt_factor_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = snowmelt_factor_textbox;
                Snowmelt_factor_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Snowmelt_factor_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                snowmelt_factor_textbox = value;
                Snowmelt_factor_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Snow_threshold_textboxRWL = new ReaderWriterLock();
        protected string snow_threshold_textbox = "";
        public string Snow_threshold_textbox
        {
            get
            {
                Snow_threshold_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = snow_threshold_textbox;
                Snow_threshold_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Snow_threshold_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                snow_threshold_textbox = value;
                Snow_threshold_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Tillfields_input_filename_textboxRWL = new ReaderWriterLock();
        protected string tillfields_input_filename_textbox = "";
        public string Tillfields_input_filename_textbox
        {
            get
            {
                Tillfields_input_filename_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = tillfields_input_filename_textbox;
                Tillfields_input_filename_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Tillfields_input_filename_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                tillfields_input_filename_textbox = value;
                Tillfields_input_filename_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Landuse_input_filename_textboxRWL = new ReaderWriterLock();
        protected string landuse_input_filename_textbox = "";
        public string Landuse_input_filename_textbox
        {
            get
            {
                Landuse_input_filename_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = landuse_input_filename_textbox;
                Landuse_input_filename_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Landuse_input_filename_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                landuse_input_filename_textbox = value;
                Landuse_input_filename_textboxRWL.ReleaseWriterLock();
            }
        }
        ReaderWriterLock Evap_input_filename_textboxRWL = new ReaderWriterLock();
        protected string evap_input_filename_textbox = "";
        public string Evap_input_filename_textbox
        {
            get
            {
                Evap_input_filename_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = evap_input_filename_textbox;
                Evap_input_filename_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Evap_input_filename_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                evap_input_filename_textbox = value;
                Evap_input_filename_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Infil_input_filename_textboxRWL = new ReaderWriterLock();
        protected string infil_input_filename_textbox = "";
        public string Infil_input_filename_textbox
        {
            get
            {
                Infil_input_filename_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = infil_input_filename_textbox;
                Infil_input_filename_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Infil_input_filename_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                infil_input_filename_textbox = value;
                Infil_input_filename_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Rain_input_filename_textboxRWL = new ReaderWriterLock();
        protected string rain_input_filename_textbox = "";
        public string Rain_input_filename_textbox
        {
            get
            {
                Rain_input_filename_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = rain_input_filename_textbox;
                Rain_input_filename_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Rain_input_filename_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                rain_input_filename_textbox = value;
                Rain_input_filename_textboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Temp_input_filename_textboxRWL = new ReaderWriterLock();
        protected string temp_input_filename_textbox = "";
        public string Temp_input_filename_textbox
        {
            get
            {
                Temp_input_filename_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = temp_input_filename_textbox;
                Temp_input_filename_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Temp_input_filename_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                temp_input_filename_textbox = value;
                Temp_input_filename_textboxRWL.ReleaseWriterLock();
            }
        }
        ReaderWriterLock TextBox_ls_transRWL = new ReaderWriterLock();
        protected string textBox_ls_trans = "";
        public string TextBox_ls_trans
        {
            get
            {
                TextBox_ls_transRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = textBox_ls_trans;
                TextBox_ls_transRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                TextBox_ls_transRWL.AcquireWriterLock(Timeout.Infinite);
                textBox_ls_trans = value;
                TextBox_ls_transRWL.ReleaseWriterLock();
            }
        }
        ReaderWriterLock TextBox_ls_cohRWL = new ReaderWriterLock();
        protected string textBox_ls_coh = "";
        public string TextBox_ls_coh
        {
            get
            {
                TextBox_ls_cohRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = textBox_ls_coh;
                TextBox_ls_cohRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                TextBox_ls_cohRWL.AcquireWriterLock(Timeout.Infinite);
                textBox_ls_coh = value;
                TextBox_ls_cohRWL.ReleaseWriterLock();
            }
        }
        ReaderWriterLock TextBox_ls_bdRWL = new ReaderWriterLock();
        protected string textBox_ls_bd = "";
        public string TextBox_ls_bd
        {
            get
            {
                TextBox_ls_bdRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = textBox_ls_bd;
                TextBox_ls_bdRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                TextBox_ls_bdRWL.AcquireWriterLock(Timeout.Infinite);
                textBox_ls_bd = value;
                TextBox_ls_bdRWL.ReleaseWriterLock();
            }
        }
        ReaderWriterLock TextBox_ls_ifrRWL = new ReaderWriterLock();
        protected string textBox_ls_ifr = "";
        public string TextBox_ls_ifr
        {
            get
            {
                TextBox_ls_ifrRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = textBox_ls_ifr;
                TextBox_ls_ifrRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                TextBox_ls_ifrRWL.AcquireWriterLock(Timeout.Infinite);
                textBox_ls_ifr = value;
                TextBox_ls_ifrRWL.ReleaseWriterLock();
            }
        }
        ReaderWriterLock Total_tillage_statuspanelRWL = new ReaderWriterLock();
        protected string total_tillage_statuspanel = "";
        public string Total_tillage_statuspanel
        {
            get
            {
                Total_tillage_statuspanelRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = total_tillage_statuspanel;
                Total_tillage_statuspanelRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Total_tillage_statuspanelRWL.AcquireWriterLock(Timeout.Infinite);
                total_tillage_statuspanel = value;
                Total_tillage_statuspanelRWL.ReleaseWriterLock();
            }
        }
        ReaderWriterLock DailyPRWL = new ReaderWriterLock();
        protected string dailyP = "";
        public string DailyP
        {
            get
            {
                DailyPRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = dailyP;
                DailyPRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                DailyPRWL.AcquireWriterLock(Timeout.Infinite);
                dailyP = value;
                DailyPRWL.ReleaseWriterLock();
            }
        }
        ReaderWriterLock DailyDRWL = new ReaderWriterLock();
        protected string dailyD = "";
        public string DailyD
        {
            get
            {
                DailyDRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = dailyD;
                DailyDRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                DailyDRWL.AcquireWriterLock(Timeout.Infinite);
                dailyD = value;
                DailyDRWL.ReleaseWriterLock();
            }
        }
        ReaderWriterLock DailyT_avgRWL = new ReaderWriterLock();
        protected string dailyT_avg = "";
        public string DailyT_avg
        {
            get
            {
                DailyT_avgRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = dailyT_avg;
                DailyT_avgRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                DailyT_avgRWL.AcquireWriterLock(Timeout.Infinite);
                dailyT_avg = value;
                DailyT_avgRWL.ReleaseWriterLock();
            }
        }
        ReaderWriterLock DailyT_minRWL = new ReaderWriterLock();
        protected string dailyT_min = "";
        public string DailyT_min
        {
            get
            {
                DailyT_minRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = dailyT_min;
                DailyT_minRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                DailyT_minRWL.AcquireWriterLock(Timeout.Infinite);
                dailyT_min = value;
                DailyT_minRWL.ReleaseWriterLock();
            }
        }
        ReaderWriterLock DailyT_maxRWL = new ReaderWriterLock();
        protected string dailyT_max = "";
        public string DailyT_max
        {
            get
            {
                DailyT_maxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = dailyT_max;
                DailyT_maxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                DailyT_maxRWL.AcquireWriterLock(Timeout.Infinite);
                dailyT_max = value;
                DailyT_maxRWL.ReleaseWriterLock();
            }
        }
        ReaderWriterLock GoogleBeginDateRWL = new ReaderWriterLock();
        protected string googleBeginDate = "";
        public string GoogleBeginDate
        {
            get
            {
                GoogleBeginDateRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = googleBeginDate;
                GoogleBeginDateRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                GoogleBeginDateRWL.AcquireWriterLock(Timeout.Infinite);
                googleBeginDate = value;
                GoogleBeginDateRWL.ReleaseWriterLock();
            }
        }
        ReaderWriterLock TextBoxAVIFileRWL = new ReaderWriterLock();
        protected string textBoxAVIFile = "";
        public string TextBoxAVIFile
        {
            get
            {
                TextBoxAVIFileRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = textBoxAVIFile;
                TextBoxAVIFileRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                TextBoxAVIFileRWL.AcquireWriterLock(Timeout.Infinite);
                textBoxAVIFile = value;
                TextBoxAVIFileRWL.ReleaseWriterLock();
            }
        }
        ReaderWriterLock Ini_CaCO3_contentRWL = new ReaderWriterLock();
        protected string ini_CaCO3_content = "";
        public string Ini_CaCO3_content
        {
            get
            {
                Ini_CaCO3_contentRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = ini_CaCO3_content;
                Ini_CaCO3_contentRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Ini_CaCO3_contentRWL.AcquireWriterLock(Timeout.Infinite);
                ini_CaCO3_content = value;
                Ini_CaCO3_contentRWL.ReleaseWriterLock();
            }
        }
        ReaderWriterLock Calibration_levels_textboxRWL = new ReaderWriterLock();
        protected string calibration_levels_textbox = "";
        public string Calibration_levels_textbox
        {
            get
            {
                Calibration_levels_textboxRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = calibration_levels_textbox;
                Calibration_levels_textboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Calibration_levels_textboxRWL.AcquireWriterLock(Timeout.Infinite);
                calibration_levels_textbox = value;
                Calibration_levels_textboxRWL.ReleaseWriterLock();
            }
        }
        ReaderWriterLock Daily_nRWL = new ReaderWriterLock();
        protected string daily_n = "";
        public string Daily_n
        {
            get
            {
                Daily_nRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = daily_n;
                Daily_nRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Daily_nRWL.AcquireWriterLock(Timeout.Infinite);
                daily_n = value;
                Daily_nRWL.ReleaseWriterLock();
            }
        }
        ReaderWriterLock Latitude_degRWL = new ReaderWriterLock();
        protected string latitude_deg = "";
        public string Latitude_deg
        {
            get
            {
                Latitude_degRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = latitude_deg;
                Latitude_degRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Latitude_degRWL.AcquireWriterLock(Timeout.Infinite);
                latitude_deg = value;
                Latitude_degRWL.ReleaseWriterLock();
            }
        }
        ReaderWriterLock Latitude_minRWL = new ReaderWriterLock();
        protected string latitude_min = "";
        public string Latitude_min
        {
            get
            {
                Latitude_minRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = latitude_min;
                Latitude_minRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Latitude_minRWL.AcquireWriterLock(Timeout.Infinite);
                latitude_min = value;
                Latitude_minRWL.ReleaseWriterLock();
            }
        }








































        //CheckBoxes

        ReaderWriterLock Water_ero_checkboxRWL = new ReaderWriterLock();
        protected bool water_ero_checkbox = false;
        public bool Water_ero_checkbox
        {
            get
            {
                Water_ero_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = water_ero_checkbox;
                Water_ero_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Water_ero_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                water_ero_checkbox = value;
                Water_ero_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Tillage_checkboxRWL = new ReaderWriterLock();
        protected bool tillage_checkbox = false;
        public bool Tillage_checkbox
        {
            get
            {
                Tillage_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = tillage_checkbox;
                Tillage_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Tillage_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                tillage_checkbox = value;
                Tillage_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Landslide_checkboxRWL = new ReaderWriterLock();
        protected bool landslide_checkbox = false;
        public bool Landslide_checkbox
        {
            get
            {
                Landslide_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = landslide_checkbox;
                Landslide_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Landslide_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                landslide_checkbox = value;
                Landslide_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Creep_active_checkboxRWL = new ReaderWriterLock();
        protected bool creep_active_checkbox = false;
        public bool Creep_active_checkbox
        {
            get
            {
                Creep_active_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = creep_active_checkbox;
                Creep_active_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Creep_active_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                creep_active_checkbox = value;
                Creep_active_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Biological_weathering_checkboxRWL = new ReaderWriterLock();
        protected bool biological_weathering_checkbox = false;
        public bool Biological_weathering_checkbox
        {
            get
            {
                Biological_weathering_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = biological_weathering_checkbox;
                Biological_weathering_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Biological_weathering_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                biological_weathering_checkbox = value;
                Biological_weathering_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Frost_weathering_checkboxRWL = new ReaderWriterLock();
        protected bool frost_weathering_checkbox = false;
        public bool Frost_weathering_checkbox
        {
            get
            {
                Frost_weathering_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = frost_weathering_checkbox;
                Frost_weathering_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Frost_weathering_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                frost_weathering_checkbox = value;
                Frost_weathering_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Tilting_active_checkboxRWL = new ReaderWriterLock();
        protected bool tilting_active_checkbox = false;
        public bool Tilting_active_checkbox
        {
            get
            {
                Tilting_active_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = tilting_active_checkbox;
                Tilting_active_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Tilting_active_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                tilting_active_checkbox = value;
                Tilting_active_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Uplift_active_checkboxRWL = new ReaderWriterLock();
        protected bool uplift_active_checkbox = false;
        public bool Uplift_active_checkbox
        {
            get
            {
                Uplift_active_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = uplift_active_checkbox;
                Uplift_active_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Uplift_active_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                uplift_active_checkbox = value;
                Uplift_active_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Soil_phys_weath_checkboxRWL = new ReaderWriterLock();
        protected bool soil_phys_weath_checkbox = false;
        public bool Soil_phys_weath_checkbox
        {
            get
            {
                Soil_phys_weath_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = soil_phys_weath_checkbox;
                Soil_phys_weath_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Soil_phys_weath_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                soil_phys_weath_checkbox = value;
                Soil_phys_weath_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Soil_chem_weath_checkboxRWL = new ReaderWriterLock();
        protected bool soil_chem_weath_checkbox = false;
        public bool Soil_chem_weath_checkbox
        {
            get
            {
                Soil_chem_weath_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = soil_chem_weath_checkbox;
                Soil_chem_weath_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Soil_chem_weath_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                soil_chem_weath_checkbox = value;
                Soil_chem_weath_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Soil_bioturb_checkboxRWL = new ReaderWriterLock();
        protected bool soil_bioturb_checkbox = false;
        public bool Soil_bioturb_checkbox
        {
            get
            {
                Soil_bioturb_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = soil_bioturb_checkbox;
                Soil_bioturb_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Soil_bioturb_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                soil_bioturb_checkbox = value;
                Soil_bioturb_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Soil_clay_transloc_checkboxRWL = new ReaderWriterLock();
        protected bool soil_clay_transloc_checkbox = false;
        public bool Soil_clay_transloc_checkbox
        {
            get
            {
                Soil_clay_transloc_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = soil_clay_transloc_checkbox;
                Soil_clay_transloc_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Soil_clay_transloc_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                soil_clay_transloc_checkbox = value;
                Soil_clay_transloc_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Soil_carbon_cycle_checkboxRWL = new ReaderWriterLock();
        protected bool soil_carbon_cycle_checkbox = false;
        public bool Soil_carbon_cycle_checkbox
        {
            get
            {
                Soil_carbon_cycle_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = soil_carbon_cycle_checkbox;
                Soil_carbon_cycle_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Soil_carbon_cycle_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                soil_carbon_cycle_checkbox = value;
                Soil_carbon_cycle_checkboxRWL.ReleaseWriterLock();
            }
        }










        ReaderWriterLock Calibration_buttonRWL = new ReaderWriterLock();
        protected bool calibration_button = false;
        public bool Calibration_button
        {
            get
            {
                Calibration_buttonRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = calibration_button;
                Calibration_buttonRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Calibration_buttonRWL.AcquireWriterLock(Timeout.Infinite);
                calibration_button = value;
                Calibration_buttonRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Sensitivity_buttonRWL = new ReaderWriterLock();
        protected bool sensitivity_button = false;
        public bool Sensitivity_button
        {
            get
            {
                Sensitivity_buttonRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = sensitivity_button;
                Sensitivity_buttonRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Sensitivity_buttonRWL.AcquireWriterLock(Timeout.Infinite);
                sensitivity_button = value;
                Sensitivity_buttonRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock UTMgridcheckboxRWL = new ReaderWriterLock();
        protected bool uTMgridcheckbox = false;
        public bool UTMgridcheckbox
        {
            get
            {
                UTMgridcheckboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = uTMgridcheckbox;
                UTMgridcheckboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                UTMgridcheckboxRWL.AcquireWriterLock(Timeout.Infinite);
                uTMgridcheckbox = value;
                UTMgridcheckboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock GoogleAnimationCheckboxRWL = new ReaderWriterLock();
        protected bool googleAnimationCheckbox = false;
        public bool GoogleAnimationCheckbox
        {
            get
            {
                GoogleAnimationCheckboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = googleAnimationCheckbox;
                GoogleAnimationCheckboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                GoogleAnimationCheckboxRWL.AcquireWriterLock(Timeout.Infinite);
                googleAnimationCheckbox = value;
                GoogleAnimationCheckboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock CheckBoxGenerateAVIFileRWL = new ReaderWriterLock();
        protected bool checkBoxGenerateAVIFile = false;
        public bool CheckBoxGenerateAVIFile
        {
            get
            {
                CheckBoxGenerateAVIFileRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = checkBoxGenerateAVIFile;
                CheckBoxGenerateAVIFileRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                CheckBoxGenerateAVIFileRWL.AcquireWriterLock(Timeout.Infinite);
                checkBoxGenerateAVIFile = value;
                CheckBoxGenerateAVIFileRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Radio_tilt_col_zeroRWL = new ReaderWriterLock();
        protected bool radio_tilt_col_zero = false;
        public bool Radio_tilt_col_zero
        {
            get
            {
                Radio_tilt_col_zeroRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = radio_tilt_col_zero;
                Radio_tilt_col_zeroRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Radio_tilt_col_zeroRWL.AcquireWriterLock(Timeout.Infinite);
                radio_tilt_col_zero = value;
                Radio_tilt_col_zeroRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Radio_tilt_row_zeroRWL = new ReaderWriterLock();
        protected bool radio_tilt_row_zero = false;
        public bool Radio_tilt_row_zero
        {
            get
            {
                Radio_tilt_row_zeroRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = radio_tilt_row_zero;
                Radio_tilt_row_zeroRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Radio_tilt_row_zeroRWL.AcquireWriterLock(Timeout.Infinite);
                radio_tilt_row_zero = value;
                Radio_tilt_row_zeroRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Radio_tilt_col_maxRWL = new ReaderWriterLock();
        protected bool radio_tilt_col_max = false;
        public bool Radio_tilt_col_max
        {
            get
            {
                Radio_tilt_col_maxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = radio_tilt_col_max;
                Radio_tilt_col_maxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Radio_tilt_col_maxRWL.AcquireWriterLock(Timeout.Infinite);
                radio_tilt_col_max = value;
                Radio_tilt_col_maxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Radio_tilt_row_maxRWL = new ReaderWriterLock();
        protected bool radio_tilt_row_max = false;
        public bool Radio_tilt_row_max
        {
            get
            {
                Radio_tilt_row_maxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = radio_tilt_row_max;
                Radio_tilt_row_maxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Radio_tilt_row_maxRWL.AcquireWriterLock(Timeout.Infinite);
                radio_tilt_row_max = value;
                Radio_tilt_row_maxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Radio_lift_row_less_thanRWL = new ReaderWriterLock();
        protected bool radio_lift_row_less_than = false;
        public bool Radio_lift_row_less_than
        {
            get
            {
                Radio_lift_row_less_thanRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = radio_lift_row_less_than;
                Radio_lift_row_less_thanRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Radio_lift_row_less_thanRWL.AcquireWriterLock(Timeout.Infinite);
                radio_lift_row_less_than = value;
                Radio_lift_row_less_thanRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Radio_lift_row_more_thanRWL = new ReaderWriterLock();
        protected bool radio_lift_row_more_than = false;
        public bool Radio_lift_row_more_than
        {
            get
            {
                Radio_lift_row_more_thanRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = radio_lift_row_more_than;
                Radio_lift_row_more_thanRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Radio_lift_row_more_thanRWL.AcquireWriterLock(Timeout.Infinite);
                radio_lift_row_more_than = value;
                Radio_lift_row_more_thanRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Radio_lift_col_less_thanRWL = new ReaderWriterLock();
        protected bool radio_lift_col_less_than = false;
        public bool Radio_lift_col_less_than
        {
            get
            {
                Radio_lift_col_less_thanRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = radio_lift_col_less_than;
                Radio_lift_col_less_thanRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Radio_lift_col_less_thanRWL.AcquireWriterLock(Timeout.Infinite);
                radio_lift_col_less_than = value;
                Radio_lift_col_less_thanRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Radio_lift_col_more_thanRWL = new ReaderWriterLock();
        protected bool radio_lift_col_more_than = false;
        public bool Radio_lift_col_more_than
        {
            get
            {
                Radio_lift_col_more_thanRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = radio_lift_col_more_than;
                Radio_lift_col_more_thanRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Radio_lift_col_more_thanRWL.AcquireWriterLock(Timeout.Infinite);
                radio_lift_col_more_than = value;
                Radio_lift_col_more_thanRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Treefall_checkboxRWL = new ReaderWriterLock();
        protected bool treefall_checkbox = false;
        public bool Treefall_checkbox
        {
            get
            {
                Treefall_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = treefall_checkbox;
                Treefall_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Treefall_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                treefall_checkbox = value;
                Treefall_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock CT_depth_decay_checkboxRWL = new ReaderWriterLock();
        protected bool cT_depth_decay_checkbox = false;
        public bool CT_depth_decay_checkbox
        {
            get
            {
                CT_depth_decay_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = cT_depth_decay_checkbox;
                CT_depth_decay_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                CT_depth_decay_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                cT_depth_decay_checkbox = value;
                CT_depth_decay_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Check_space_soildepthRWL = new ReaderWriterLock();
        protected bool check_space_soildepth = false;
        public bool Check_space_soildepth
        {
            get
            {
                Check_space_soildepthRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = check_space_soildepth;
                Check_space_soildepthRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Check_space_soildepthRWL.AcquireWriterLock(Timeout.Infinite);
                check_space_soildepth = value;
                Check_space_soildepthRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Check_space_landuseRWL = new ReaderWriterLock();
        protected bool check_space_landuse = false;
        public bool Check_space_landuse
        {
            get
            {
                Check_space_landuseRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = check_space_landuse;
                Check_space_landuseRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Check_space_landuseRWL.AcquireWriterLock(Timeout.Infinite);
                check_space_landuse = value;
                Check_space_landuseRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Check_time_landuseRWL = new ReaderWriterLock();
        protected bool check_time_landuse = false;
        public bool Check_time_landuse
        {
            get
            {
                Check_time_landuseRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = check_time_landuse;
                Check_time_landuseRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Check_time_landuseRWL.AcquireWriterLock(Timeout.Infinite);
                check_time_landuse = value;
                Check_time_landuseRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Check_space_evapRWL = new ReaderWriterLock();
        protected bool check_space_evap = false;
        public bool Check_space_evap
        {
            get
            {
                Check_space_evapRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = check_space_evap;
                Check_space_evapRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Check_space_evapRWL.AcquireWriterLock(Timeout.Infinite);
                check_space_evap = value;
                Check_space_evapRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Check_time_evapRWL = new ReaderWriterLock();
        protected bool check_time_evap = false;
        public bool Check_time_evap
        {
            get
            {
                Check_time_evapRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = check_time_evap;
                Check_time_evapRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Check_time_evapRWL.AcquireWriterLock(Timeout.Infinite);
                check_time_evap = value;
                Check_time_evapRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Check_space_infilRWL = new ReaderWriterLock();
        protected bool check_space_infil = false;
        public bool Check_space_infil
        {
            get
            {
                Check_space_infilRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = check_space_infil;
                Check_space_infilRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Check_space_infilRWL.AcquireWriterLock(Timeout.Infinite);
                check_space_infil = value;
                Check_space_infilRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Check_time_infilRWL = new ReaderWriterLock();
        protected bool check_time_infil = false;
        public bool Check_time_infil
        {
            get
            {
                Check_time_infilRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = check_time_infil;
                Check_time_infilRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Check_time_infilRWL.AcquireWriterLock(Timeout.Infinite);
                check_time_infil = value;
                Check_time_infilRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Check_space_rainRWL = new ReaderWriterLock();
        protected bool check_space_rain = false;
        public bool Check_space_rain
        {
            get
            {
                Check_space_rainRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = check_space_rain;
                Check_space_rainRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Check_space_rainRWL.AcquireWriterLock(Timeout.Infinite);
                check_space_rain = value;
                Check_space_rainRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Check_time_rainRWL = new ReaderWriterLock();
        protected bool check_time_rain = false;
        public bool Check_time_rain
        {
            get
            {
                Check_time_rainRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = check_time_rain;
                Check_time_rainRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Check_time_rainRWL.AcquireWriterLock(Timeout.Infinite);
                check_time_rain = value;
                Check_time_rainRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Check_time_TRWL = new ReaderWriterLock();
        protected bool check_time_T = false;
        public bool Check_time_T
        {
            get
            {
                Check_time_TRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = check_time_T;
                Check_time_TRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Check_time_TRWL.AcquireWriterLock(Timeout.Infinite);
                check_time_T = value;
                Check_time_TRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Only_waterflow_checkboxRWL = new ReaderWriterLock();
        protected bool only_waterflow_checkbox = false;
        public bool Only_waterflow_checkbox
        {
            get
            {
                Only_waterflow_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = only_waterflow_checkbox;
                Only_waterflow_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Only_waterflow_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                only_waterflow_checkbox = value;
                Only_waterflow_checkboxRWL.ReleaseWriterLock();
            }
        }


        ReaderWriterLock Version_lux_checkboxRWL = new ReaderWriterLock();
        protected bool version_lux_checkbox = false;
        public bool Version_lux_checkbox
        {
            get
            {
                Version_lux_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = version_lux_checkbox;
                Version_lux_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Version_lux_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                version_lux_checkbox = value;
                Version_lux_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Solifluction_checkboxRWL = new ReaderWriterLock();
        protected bool solifluction_checkbox = false;
        public bool Solifluction_checkbox
        {
            get
            {
                Solifluction_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = solifluction_checkbox;
                Solifluction_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Solifluction_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                solifluction_checkbox = value;
                Solifluction_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Decalcification_checkboxRWL = new ReaderWriterLock();
        protected bool decalcification_checkbox = false;
        public bool Decalcification_checkbox
        {
            get
            {
                Decalcification_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = decalcification_checkbox;
                Decalcification_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Decalcification_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                decalcification_checkbox = value;
                Decalcification_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Daily_waterRWL = new ReaderWriterLock();
        protected bool daily_water = false;
        public bool Daily_water
        {
            get
            {
                Daily_waterRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = daily_water;
                Daily_waterRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Daily_waterRWL.AcquireWriterLock(Timeout.Infinite);
                daily_water = value;
                Daily_waterRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Ct_JagercikovaRWL = new ReaderWriterLock();
        protected bool ct_Jagercikova = false;
        public bool Ct_Jagercikova
        {
            get
            {
                Ct_JagercikovaRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = ct_Jagercikova;
                Ct_JagercikovaRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Ct_JagercikovaRWL.AcquireWriterLock(Timeout.Infinite);
                ct_Jagercikova = value;
                Ct_JagercikovaRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock View_maps_checkboxRWL = new ReaderWriterLock();
        protected bool view_maps_checkbox = false;
        public bool View_maps_checkbox
        {
            get
            {
                View_maps_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = view_maps_checkbox;
                View_maps_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                View_maps_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                view_maps_checkbox = value;
                View_maps_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Final_output_checkboxRWL = new ReaderWriterLock();
        protected bool final_output_checkbox = false;
        public bool Final_output_checkbox
        {
            get
            {
                Final_output_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = final_output_checkbox;
                Final_output_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Final_output_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                final_output_checkbox = value;
                Final_output_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Regular_output_checkboxRWL = new ReaderWriterLock();
        protected bool regular_output_checkbox = false;
        public bool Regular_output_checkbox
        {
            get
            {
                Regular_output_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = regular_output_checkbox;
                Regular_output_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Regular_output_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                regular_output_checkbox = value;
                Regular_output_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Altitude_output_checkboxRWL = new ReaderWriterLock();
        protected bool altitude_output_checkbox = false;
        public bool Altitude_output_checkbox
        {
            get
            {
                Altitude_output_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = altitude_output_checkbox;
                Altitude_output_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Altitude_output_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                altitude_output_checkbox = value;
                Altitude_output_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Creep_CheckboxRWL = new ReaderWriterLock();
        protected bool creep_Checkbox = false;
        public bool Creep_Checkbox
        {
            get
            {
                Creep_CheckboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = creep_Checkbox;
                Creep_CheckboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Creep_CheckboxRWL.AcquireWriterLock(Timeout.Infinite);
                creep_Checkbox = value;
                Creep_CheckboxRWL.ReleaseWriterLock();
            }
        }


        ReaderWriterLock Check_space_till_fieldsRWL = new ReaderWriterLock();
        protected bool check_space_till_fields = false;
        public bool Check_space_till_fields
        {
            get
            {
                Check_space_till_fieldsRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = check_space_till_fields;
                Check_space_till_fieldsRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Check_space_till_fieldsRWL.AcquireWriterLock(Timeout.Infinite);
                check_space_till_fields = value;
                Check_space_till_fieldsRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Fill_sinks_before_checkboxRWL = new ReaderWriterLock();
        protected bool fill_sinks_before_checkbox = false;
        public bool Fill_sinks_before_checkbox
        {
            get
            {
                Fill_sinks_before_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = fill_sinks_before_checkbox;
                Fill_sinks_before_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Fill_sinks_before_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                fill_sinks_before_checkbox = value;
                Fill_sinks_before_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Creep_testingRWL = new ReaderWriterLock();
        protected bool creep_testing = false;
        public bool Creep_testing
        {
            get
            {
                Creep_testingRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = creep_testing;
                Creep_testingRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Creep_testingRWL.AcquireWriterLock(Timeout.Infinite);
                creep_testing = value;
                Creep_testingRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Soildepth_output_checkboxRWL = new ReaderWriterLock();
        protected bool soildepth_output_checkbox = false;
        public bool Soildepth_output_checkbox
        {
            get
            {
                Soildepth_output_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = soildepth_output_checkbox;
                Soildepth_output_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Soildepth_output_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                soildepth_output_checkbox = value;
                Soildepth_output_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Alt_change_output_checkboxRWL = new ReaderWriterLock();
        protected bool alt_change_output_checkbox = false;
        public bool Alt_change_output_checkbox
        {
            get
            {
                Alt_change_output_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = alt_change_output_checkbox;
                Alt_change_output_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Alt_change_output_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                alt_change_output_checkbox = value;
                Alt_change_output_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Water_output_checkboxRWL = new ReaderWriterLock();
        protected bool water_output_checkbox = false;
        public bool Water_output_checkbox
        {
            get
            {
                Water_output_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = water_output_checkbox;
                Water_output_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Water_output_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                water_output_checkbox = value;
                Water_output_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Depressions_output_checkboxRWL = new ReaderWriterLock();
        protected bool depressions_output_checkbox = false;
        public bool Depressions_output_checkbox
        {
            get
            {
                Depressions_output_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = depressions_output_checkbox;
                Depressions_output_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Depressions_output_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                depressions_output_checkbox = value;
                Depressions_output_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Diagnostic_output_checkboxRWL = new ReaderWriterLock();
        protected bool diagnostic_output_checkbox = false;
        public bool Diagnostic_output_checkbox
        {
            get
            {
                Diagnostic_output_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = diagnostic_output_checkbox;
                Diagnostic_output_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Diagnostic_output_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                diagnostic_output_checkbox = value;
                Diagnostic_output_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock All_process_output_checkboxRWL = new ReaderWriterLock();
        protected bool all_process_output_checkbox = false;
        public bool All_process_output_checkbox
        {
            get
            {
                All_process_output_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = all_process_output_checkbox;
                All_process_output_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                All_process_output_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                all_process_output_checkbox = value;
                All_process_output_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Annual_output_checkboxRWL = new ReaderWriterLock();
        protected bool annual_output_checkbox = false;
        public bool Annual_output_checkbox
        {
            get
            {
                Annual_output_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = annual_output_checkbox;
                Annual_output_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Annual_output_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                annual_output_checkbox = value;
                Annual_output_checkboxRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Fill_sinks_during_checkboxRWL = new ReaderWriterLock();
        protected bool fill_sinks_during_checkbox = false;
        public bool Fill_sinks_during_checkbox
        {
            get
            {
                Fill_sinks_during_checkboxRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = fill_sinks_during_checkbox;
                Fill_sinks_during_checkboxRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Fill_sinks_during_checkboxRWL.AcquireWriterLock(Timeout.Infinite);
                fill_sinks_during_checkbox = value;
                Fill_sinks_during_checkboxRWL.ReleaseWriterLock();
            }
        }







        ReaderWriterLock UTMsouthcheckRWL = new ReaderWriterLock();
        protected bool uTMsouthcheck = false;
        public bool UTMsouthcheck
        {
            get
            {
                UTMsouthcheckRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = uTMsouthcheck;
                UTMsouthcheckRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                UTMsouthcheckRWL.AcquireWriterLock(Timeout.Infinite);
                uTMsouthcheck = value;
                UTMsouthcheckRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Check_time_till_fieldsRWL = new ReaderWriterLock();
        protected bool check_time_till_fields = false;
        public bool Check_time_till_fields
        {
            get
            {
                Check_time_till_fieldsRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = check_time_till_fields;
                Check_time_till_fieldsRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Check_time_till_fieldsRWL.AcquireWriterLock(Timeout.Infinite);
                check_time_till_fields = value;
                Check_time_till_fieldsRWL.ReleaseWriterLock();
            }
        }

        ReaderWriterLock Check_scaling_daily_weatherRWL = new ReaderWriterLock();
        protected bool check_scaling_daily_weather = false;
        public bool Check_scaling_daily_weather
        {
            get
            {
                Check_scaling_daily_weatherRWL.AcquireReaderLock(Timeout.Infinite);
                bool temp = check_scaling_daily_weather;
                Check_scaling_daily_weatherRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                Check_scaling_daily_weatherRWL.AcquireWriterLock(Timeout.Infinite);
                check_scaling_daily_weather = value;
                Check_scaling_daily_weatherRWL.ReleaseWriterLock();
            }
        }



        
    }
}
