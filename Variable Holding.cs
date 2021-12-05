using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;

namespace LORICAVariables
{
    public class GUIVariables
    {
        public GUIVariables(DelUpdateAllFields UAF, DelUpdateStatusPannel UPS, DelUpdateTimePannel UTP, DelUpdateVariables UF)
        {
            updateAllFields = UAF;
            updateStatusPannel = UPS;
            updateTimePannel = UTP;
            updateVariables = UF;

            updateVariables();
        }

        public delegate void DelUpdateAllFields(); //Pushes Variables to GUI
        public delegate void DelUpdateStatusPannel(); //Updates GUI Status Pannel
        public delegate void DelUpdateTimePannel(); //Updates GUI Status Pannel
        public delegate void DelUpdateVariables(); //Pulls Values from GUI

        DelUpdateAllFields updateAllFields;
        DelUpdateStatusPannel updateStatusPannel;
        DelUpdateTimePannel updateTimePannel;
        DelUpdateVariables updateVariables;


        //public delegate double Delbulk_Density_Calc(double, double , double, double, double, double, double, double); //Pulls Values from GUI
        //public Delbulk_Density_Calc bulk_density_calc;

        //if none of these fields get additional instructions, remove and rename the delegates
        //Ex: updateAllFields becomes UpdateAllFields
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




    }
}
