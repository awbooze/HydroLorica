

// HydroLorica v1.0 by W.M. van der Meij and A.J.A.M. Temme, 2020
// Based on Lorica by Vanwalleghem and Temme 2016
// Based on MILESD 2011 by Vanwalleghem et al and on LAPSUS by Temme, Schoorl and colleagues (2006-2011)
// 
// Credits to T.J. Coulthard for interface coding template (the CAESAR model, www.coulthard.org.uk)

//This program is free software; you can redistribute it and/or modify it under the terms of the 
//GNU General Public License as published by the Free Software Foundation;  
//This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
//without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
//See the GNU General Public License (http://www.gnu.org/copyleft/gpl.html) for more details. 
//You should have received a copy of the GNU General Public License along with this program; 
//if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, 
//MA 02110-1301, USA.

// June 2020


using System;
using System.Collections;
using System.Collections.Generic;
using System.Collections.Concurrent;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.IO;
using System.Net;
using System.Text;
using System.Windows.Forms;
using System.Xml;
using System.Runtime.InteropServices;
using System.Diagnostics;
using System.Numerics;
using MathNet.Numerics;
using MathNet.Numerics.IntegralTransforms;

using LORICAVariables;

namespace LORICA4
{
    /// <summary>
    /// Main LORICA interface
    /// </summary>
    public class Mother_form : System.Windows.Forms.Form
    {
        [DllImport("msvcrt")]
        static extern int _getch();

        private Task StartThread; //Added by Parallelization Team
        private GUIVariables guiVariables;
        Simulation MainSimulation = null;

        #region global interface parameters
        private TabControl Process_tabs;
        private TabPage Water;
        private TextBox parameter_n_textbox;
        private TextBox parameter_conv_textbox;
        private TextBox parameter_K_textbox;
        private TextBox parameter_m_textbox;
        private CheckBox only_waterflow_checkbox;
        private PictureBox pictureBox1;
        private Label label12;
        private Label label11;
        private Label label10;
        private Label label9;
        private CheckBox Water_ero_checkbox;
        private TabPage Tillage;
        private PictureBox pictureBox2;
        private Label label20;
        private Label trte;
        private TextBox parameter_tillage_constant_textbox;
        private TextBox parameter_ploughing_depth_textbox;
        private CheckBox Tillage_checkbox;
        private TabPage Creeper;
        private PictureBox pictureBox3;
        private Label label19;
        private TextBox parameter_diffusivity_textbox;
        private CheckBox creep_active_checkbox;
        private Label label36;
        private RadioButton radio_ls_fraction;
        private RadioButton radio_ls_absolute;
        private Label label35;
        private Label label34;
        private TextBox text_ls_rel_rain_intens;
        private TextBox textBox_ls_trans;
        private TextBox textBox_ls_bd;
        private TextBox textBox_ls_ifr;
        private TextBox textBox_ls_coh;
        private TextBox text_ls_abs_rain_intens;
        private Label label32;
        private Label label31;
        private Label label30;
        private Label label22;
        private Label label18;
        private PictureBox pictureBox4;
        private CheckBox Landslide_checkbox;
        private TabPage Solifluction;
        private PictureBox pictureBox5;
        private CheckBox Solifluction_checkbox;
        private TabPage Rock_weathering;
        private PictureBox pictureBox6;
        private GroupBox groupBox10;
        private CheckBox Frost_weathering_checkbox;
        private GroupBox groupBox9;
        private TextBox parameter_k1_textbox;
        private Label label24;
        private Label label26;
        private Label label27;
        private Label label28;
        private TextBox parameter_k2_textbox;
        private TextBox parameter_Pa_textbox;
        private TextBox parameter_P0_textbox;
        private Label label21;
        private CheckBox Biological_weathering_checkbox;
        private TabPage Tectonics;
        private GroupBox groupBox14;
        private GroupBox groupBox16;
        private TextBox text_lift_col_less;
        private TextBox text_lift_col_more;
        private TextBox text_lift_row_less;
        private TextBox text_lift_row_more;
        private RadioButton radio_lift_col_less_than;
        private RadioButton radio_lift_row_more_than;
        private RadioButton radio_lift_col_more_than;
        private RadioButton radio_lift_row_less_than;
        private TextBox Uplift_rate_textbox;
        private CheckBox uplift_active_checkbox;
        private Label label39;
        private GroupBox groupBox4;
        private Label label38;
        private TextBox Tilting_rate_textbox;
        private GroupBox groupBox15;
        private RadioButton radio_tilt_col_max;
        private RadioButton radio_tilt_row_zero;
        private RadioButton radio_tilt_col_zero;
        private RadioButton radio_tilt_row_max;
        private CheckBox tilting_active_checkbox;
        private TabPage tabPage1;
        private TabControl tabControl2;
        private TabPage physical;
        private TabPage chemical;
        private TabPage clay;
        private TabPage bioturbation;
        private CheckBox soil_phys_weath_checkbox;
        private CheckBox soil_chem_weath_checkbox;
        private CheckBox soil_clay_transloc_checkbox;
        private CheckBox soil_bioturb_checkbox;
        private TextBox upper_particle_fine_clay_textbox;
        private TextBox upper_particle_clay_textbox;
        private TextBox upper_particle_silt_textbox;
        private TextBox upper_particle_sand_textbox;
        private TextBox upper_particle_coarse_textbox;
        private TextBox physical_weath_constant2;
        private TextBox physical_weath_constant1;
        private TextBox Physical_weath_C1_textbox;
        private TextBox chem_weath_specific_coefficient_textbox;
        private TextBox chem_weath_depth_constant_textbox;
        private TextBox chem_weath_rate_constant_textbox;
        private TextBox specific_area_fine_clay_textbox;
        private TextBox specific_area_clay_textbox;
        private TextBox specific_area_silt_textbox;
        private TextBox specific_area_sand_textbox;
        private TextBox specific_area_coarse_textbox;
        private TextBox clay_neoform_C2_textbox;
        private TextBox clay_neoform_C1_textbox;
        private TextBox clay_neoform_constant_textbox;
        private TextBox eluviation_coefficient_textbox;
        private TextBox maximum_eluviation_textbox;
        private TextBox bioturbation_depth_decay_textbox;
        private TextBox potential_bioturbation_textbox;
        private TabPage carbon;
        private TextBox carbon_y_depth_decay_textbox;
        private TextBox carbon_humification_fraction_textbox;
        private TextBox carbon_depth_decay_textbox;
        private TextBox carbon_input_textbox;
        private CheckBox soil_carbon_cycle_checkbox;
        private TextBox carbon_o_twi_decay_textbox;
        private TextBox carbon_y_twi_decay_textbox;
        private TextBox carbon_o_depth_decay_textbox;
        private TextBox carbon_o_decomp_rate_textbox;
        private TextBox carbon_y_decomp_rate_textbox;


        private System.Windows.Forms.MainMenu mainMenu1;
        private System.Windows.Forms.MenuItem menuItemConfigFile;
        private System.Windows.Forms.MenuItem menuItemConfigFileOpen;
        private System.Windows.Forms.MenuItem menuItemConfigFileSave;
        private System.Windows.Forms.MenuItem menuItemConfigFileSaveAs;
        private System.Windows.Forms.StatusBar statusBar1;
        private System.Windows.Forms.StatusBarPanel TimeStatusPanel;
        private System.Windows.Forms.StatusBarPanel ProcessStatusPanel;
        private System.Windows.Forms.StatusBarPanel InfoStatusPanel;
        private System.Windows.Forms.Button start_button;
        private System.Windows.Forms.Button End_button;
        private System.Windows.Forms.ToolTip toolTip1;
        private Smallwisdom.Windows.Forms.ZoomPanImageBox Mapwindow;
        private MenuItem Mapselector;
        private MenuItem Menu_map_total_sediment;
        private MenuItem Menu_map_waterflow;
        private TrackBar trackBar1;
        private Label label61;
        private GroupBox map_controls;
        private ComboBox comboBox1;
        private Label label62;
        private Label label63;
        private TrackBar trackBar2;
        private CheckBox View_tabs_checkbox;
        private OpenFileDialog openFileDialog1;
        private Label label1;
        private Button button6;
        private TextBox textBox1;
        private TextBox textBox2;
        private Label label2;
        private StatusBarPanel out_sed_statuspanel;
        private StatusBarPanel total_tillage_statuspanel;
        private TabPage Output;
        private GroupBox groupBox6;
        private GroupBox groupBox1;
        private CheckBox water_output_checkbox;
        private CheckBox depressions_output_checkbox;
        private CheckBox all_process_output_checkbox;
        private CheckBox Soildepth_output_checkbox;
        private CheckBox Alt_change_output_checkbox;
        private CheckBox Altitude_output_checkbox;
        private CheckedListBox checkedListBox1;
        private CheckBox Regular_output_checkbox;
        private CheckBox Final_output_checkbox;
        private TextBox Box_years_output;
        private GroupBox groupBox5;
        private GroupBox UTMgroupBox;
        private TextBox textBox6;
        private CheckBox UTMsouthcheck;
        private TextBox UTMzonebox;
        private CheckBox UTMgridcheckbox;
        private TextBox textBox4;
        private TextBox googleBeginDate;
        private TextBox googAnimationSaveInterval;
        private TextBox googleAnimationTextBox;
        private TextBox saveintervalbox;
        private TextBox textBoxAVIFile;
        private Label label78;
        private Label label79;
        private CheckBox googleAnimationCheckbox;
        private Label label33;
        private CheckBox checkBoxGenerateAVIFile;
        private TabPage Run;
        private GroupBox groupBox7;
        private RadioButton runs_checkbox;
        private Label label16;
        private TextBox Number_runs_textbox;
        private TabPage Input;
        private TextBox tillfields_constant_textbox;
        private TextBox tillfields_input_filename_textbox;
        private TextBox evap_constant_value_box;
        private TextBox evap_input_filename_textbox;
        private TextBox infil_constant_value_box;
        private TextBox infil_input_filename_textbox;
        private TextBox rainfall_constant_value_box;
        private TextBox landuse_constant_value_box;
        private TextBox soildepth_constant_value_box;
        private TextBox landuse_input_filename_textbox;
        private TextBox soildepth_input_filename_textbox;
        private TextBox rain_input_filename_textbox;
        private TextBox dtm_input_filename_textbox;
        private GroupBox groupBox8;
        private CheckBox fill_sinks_before_checkbox;
        private CheckBox check_space_evap;
        private CheckBox check_space_infil;
        private CheckBox check_space_rain;
        private CheckBox check_space_till_fields;
        private CheckBox check_space_landuse;
        private CheckBox check_space_soildepth;
        private Label label17;
        private Label label15;
        private Label label14;
        private Label label7;
        private Label label5;
        private Label label4;
        private Label label3;
        private Label label25;
        private Label label23;
        private TabPage Processes;
        private Button graphicToGoogleEarthButton;
        private CheckBox Creep_Checkbox;
        private TabControl tabControl1;
        private GroupBox groupBox12;
        private GroupBox groupBox11;
        private Label label8;
        private RadioButton annual_output_checkbox;
        private RadioButton cumulative_output_checkbox;
        private GroupBox groupBox13;
        private CheckBox fill_sinks_during_checkbox;
        private CheckBox check_space_DTM;
        private CheckBox check_time_evap;
        private CheckBox check_time_infil;
        private CheckBox check_time_rain;
        private CheckBox check_time_till_fields;
        private CheckBox check_time_landuse;
        private Label label29;
        private Button explain_input_button;
        private MenuItem Menu_About_box;
        private MenuItem Menu_map_tillage;
        private MenuItem Menu_map_water_ero;
        private MenuItem Menu_map_creep;
        private MenuItem Menu_map_weathering;
        private Button timeseries_form_button;
        private Button button2;
        private Button button3;
        private Button button1;
        private CheckBox view_maps_checkbox;
        private MenuItem Menu_map_landsliding;
        private Label label37;
        private TextBox outputcode_textbox;
        private CheckBox diagnostic_output_checkbox;
        private GroupBox groupBox3;
        private Button landuse_determinator_button;
        private MenuItem Menu_map_critical_rainfall;
        #endregion

        #region global model parameters
        //AviWriter aw; // <JMW 20041018>
        private Label label87;
        private TextBox selectivity_constant_textbox;
        private TextBox bio_protection_constant_textbox;
        private TextBox erosion_threshold_textbox;
        private TextBox rock_protection_constant_textbox;
        private Label label90;
        private Label label91;
        private Label label92;
        private Label label88;
        private Button soil_specify_button;
        private CheckBox Ik_ben_Marijn;
        private CheckBox CT_depth_decay_checkbox;
        private TextBox ct_depth_decay;
        private CheckBox calibration;
        private CheckBox creep_testing;
        private ComboBox rockweath_method;
        private CheckBox daily_water;
        private TabPage decalcification;
        private CheckBox decalcification_checkbox;
        private Label label94;
        private TextBox ini_CaCO3_content;
        private TabPage treefall;
        private CheckBox treefall_checkbox;
        private Label label98;
        private TextBox temp_input_filename_textbox;
        private TextBox temp_constant_value_box;
        private Label label99;
        private CheckBox check_time_T;
        private TabPage tabPage2;
        private Label label105;
        private TextBox snowmelt_factor_textbox;
        private Label label104;
        private TextBox latitude_min;
        private Label label103;
        private TextBox latitude_deg;
        private Label label100;
        private Label label101;
        private Label label102;
        private TextBox dailyT_min;
        private TextBox dailyT_max;
        private TextBox dailyT_avg;
        private Label label97;
        private TextBox daily_n;
        private Label label96;
        private Label label93;
        private Label label89;
        private Label label40;
        private TextBox dailyET0;
        private TextBox dailyD;
        private TextBox dailyP;
        private Label label106;
        private TextBox snow_threshold_textbox;
        private TextBox ct_v0_Jagercikova;
        private TextBox ct_dd_Jagercikova;
        private System.Windows.Forms.Timer timer1;
        private Label label109;
        private Label label108;
        private CheckBox ct_Jagercikova;
        private CheckBox check_scaling_daily_weather;
        private TextBox tf_D;
        private Label label95;
        private Label label107;
        private TextBox tf_W;
        private TextBox tf_growth;
        private Label label110;
        private TextBox tf_age;
        private Label label111;
        private TextBox tf_freq;
        private Label label112;
        private GroupBox groupBox2;
        private Label label118;
        private Label label117;
        private Label label115;
        private Label label114;
        private RadioButton Sensitivity_button;
        private RadioButton Calibration_button;
        private Label label113;
        private TextBox calibration_ratios_textbox;
        private TextBox calibration_levels_textbox;
        private Label label116;
        private TextBox calibration_ratio_reduction_parameter_textbox;
        private Label label119;
        private Label label120;
        private CheckBox version_lux_checkbox;
        private Button button4;
        public static int plotType = 0;
        public static double magnifyValue = 0;
        public static int updateClick = 0;
        private double[] zoomFactor = { .25, .33, .50, .66, .80, 1, 1.25, 1.5, 2.0, 2.5, 3.0 };
        private double[] contrastFactor = { 1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3 };
        private double contrastMultiplier = 1;


        private System.ComponentModel.IContainer components;

        double[,,]     //3D matrices for properties of soil layers in different x y (x,y,z)
                    layerthickness_m,         // : thickness in m 
                    young_SOM_kg,         // : OM mass in kgrams (per voxel = layer * thickness)
                    old_SOM_kg,         // : OM mass in kgrams (per voxel = layer * thickness) 
                    bulkdensity;            // : bulkdensity in kg/m3 (over the voxel = layer * thickness)


        



        int calibration_length; // to assess if a process has to be calibrated mvdm, >1 is more runs

        //sinks and depression parameters:
        //the constant values below may have to be increased for large or strange landscapes and studies
        const int numberofsinks = 10000;           // run the program once to find out the number of sinks. The exact number and any higher number will do....
        
        const double root = 7.07;

        int     i,
                j,
                ii,
                jj,
                z,
                graphics_scale = 2; /*
                
                twoequals,      // the -equals are counters for different types of sinks
                threeequals,
                moreequals,
                round,
                s_ch,
                numtel,
                S1_error,
                S2_error,
                cell_lock,
                temp,
                num_str,
                matrixresult,
                er_ifile,
                flat,
                depressionreallyready,
                rememberrow,
                remembercol,
                search,
                twice_dtm_fill,
                three_dtm_fill,
                num_out,
                //ntr,				//WVG 22-10-2010 number of rows (timesteps) in profile timeseries matrices			//unused
                cross1, 			//WVG 22-10-2010 rows (or in the future columns) of which profiles are wanted
                cross2,
                cross3,
                GlobalMethods.nr,
                GlobalMethods.nc;*/

        private void rain_input_filename_textbox_TextChanged_1(object sender, EventArgs e)
        {

        }

        private void tillfields_input_filename_textbox_TextChanged_1(object sender, EventArgs e)
        {

        }

        private void label93_Click(object sender, EventArgs e)
        {

        }

        private void comboBox2_SelectedIndexChanged(object sender, EventArgs e)
        {

        }


        private void soil_chem_weath_checkbox_CheckedChanged(object sender, EventArgs e)
        {

        }

        private void decalcification_checkbox_CheckedChanged(object sender, EventArgs e)
        {

        }

        private void checkBox1_CheckedChanged_2(object sender, EventArgs e)
        {

        }

        private void label95_Click(object sender, EventArgs e)
        {

        }

        private void label11_Click(object sender, EventArgs e)
        {

        }

        private void label10_Click(object sender, EventArgs e)
        {

        }

        private void button4_Click(object sender, EventArgs e)
        {
            Debug.Write(" merely_calculating_derivatives");
            GlobalMethods.merely_calculating_derivatives = true;
            try { GlobalMethods.calculate_terrain_derivatives(); MessageBox.Show("terrain derivatives calculation succeeded"); }
            catch { MessageBox.Show("terrain derivatives calculation failed"); }
        }


        private void radioButton1_CheckedChanged(object sender, EventArgs e)
        {
            if (Calibration_button.Checked == true) { Sensitivity_button.Checked = false; }
        }

        private void radioButton2_CheckedChanged(object sender, EventArgs e)
        {
            if (Sensitivity_button.Checked == true) { Calibration_button.Checked = false; }
        }

        private void dailyT_max_TextChanged(object sender, EventArgs e)
        {
            OpenFileDialog openFileDialog1 = new OpenFileDialog();

            openFileDialog1.InitialDirectory = GlobalMethods.Workdir;
            openFileDialog1.FilterIndex = 1;
            openFileDialog1.RestoreDirectory = false;

            dailyT_max.Text = GetDialogFileName(openFileDialog1);

        }

        private void dailyT_min_TextChanged(object sender, EventArgs e)
        {
            OpenFileDialog openFileDialog1 = new OpenFileDialog();

            openFileDialog1.InitialDirectory = GlobalMethods.Workdir;
            openFileDialog1.FilterIndex = 1;
            openFileDialog1.RestoreDirectory = false;

            dailyT_min.Text = GetDialogFileName(openFileDialog1);

        }

        private void dailyT_avg_TextChanged(object sender, EventArgs e)
        {
            OpenFileDialog openFileDialog1 = new OpenFileDialog();

            openFileDialog1.InitialDirectory = GlobalMethods.Workdir;
            openFileDialog1.FilterIndex = 1;
            openFileDialog1.RestoreDirectory = false;

            dailyT_avg.Text = GetDialogFileName(openFileDialog1);

        }

        private void dailyP_TextChanged_1(object sender, EventArgs e)
        {

        }

        private void textBox3_TextChanged_2(object sender, EventArgs e)
        {

        }

        private void label103_Click(object sender, EventArgs e)
        {

        }

        private void textBox3_TextChanged_1(object sender, EventArgs e)
        {

        }

        private void label98_Click(object sender, EventArgs e)
        {

        }

        private void daily_water_CheckedChanged(object sender, EventArgs e)
        {
            dailyP.Enabled = (daily_water.CheckState == CheckState.Checked);
            dailyET0.Enabled = (daily_water.CheckState == CheckState.Checked);
            dailyD.Enabled = (daily_water.CheckState == CheckState.Checked);
            daily_n.Enabled = (daily_water.CheckState == CheckState.Checked);
            dailyT_avg.Enabled = (daily_water.CheckState == CheckState.Checked);
            dailyT_min.Enabled = (daily_water.CheckState == CheckState.Checked);
            dailyT_max.Enabled = (daily_water.CheckState == CheckState.Checked);
            temp_constant_value_box.Enabled = (daily_water.CheckState == CheckState.Checked);
            temp_input_filename_textbox.Enabled = (daily_water.CheckState == CheckState.Checked);
            check_time_T.Enabled = (daily_water.CheckState == CheckState.Checked);
            latitude_deg.Enabled = (daily_water.CheckState == CheckState.Checked);
            latitude_min.Enabled = (daily_water.CheckState == CheckState.Checked);
            snowmelt_factor_textbox.Enabled = (daily_water.CheckState == CheckState.Checked);
            snow_threshold_textbox.Enabled = (daily_water.CheckState == CheckState.Checked);

        }

        private void label97_Click(object sender, EventArgs e)
        {

        }

        private void label96_Click(object sender, EventArgs e)
        {

        }


        

        //  variables for displaying purposes // straight from Tom Coulthard
        double hue = 360.0;		// Ranges between 0 and 360 degrees
        double sat = 0.90;		// Ranges between 0 and 1.0 (where 1 is 100%)
        double val = 1.0;		// Ranges between 0 and 1.0 (where 1 is 100%)
        double red = 0.0;
        double green = 0.0;
        double blue = 0.0;

        string basetext = "LORICA Landscape Evolution Model";
        string cfgname = null;  //Config file name
        // string GlobalMethods.Workdir = "D:\\PhD\\projects\\1g_basic LORICA development\\";
        #endregion

        public Mother_form()
        {
            //
            // Required for Windows Form Designer support
            //
            InitializeComponent();

            //
            // TODO: Add any constructor code after InitializeComponent call
            //
        }

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        protected override void Dispose(bool disposing)
        {
            if (disposing)
            {
                if (components != null)
                {
                    components.Dispose();
                }
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code
        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            guiVariables = new GUIVariables(UpdateAllFields, UpdateStatusPannel, UpdateTimePannel, DelUpdateVariables, UpdateTimeSeries, UpdateProfile, UpdateLanduse_determinator, UpdateSoildata, draw_map);


            GlobalMethods.Workdir = "D:\\PhD\\projects\\2g_clorpt effects on soil landscape diversity\\";

            this.components = new System.ComponentModel.Container();
            System.Windows.Forms.Label label6;
            System.Windows.Forms.TabPage Landsliding;
            System.ComponentModel.ComponentResourceManager resources = new System.ComponentModel.ComponentResourceManager(typeof(Mother_form));
            System.Windows.Forms.Label label41;
            System.Windows.Forms.Label label42;
            System.Windows.Forms.Label label43;
            System.Windows.Forms.Label label44;
            System.Windows.Forms.Label label45;
            System.Windows.Forms.Label label46;
            System.Windows.Forms.Label label47;
            System.Windows.Forms.Label label48;
            System.Windows.Forms.Label label49;
            System.Windows.Forms.Label label50;
            System.Windows.Forms.Label label51;
            System.Windows.Forms.Label label52;
            System.Windows.Forms.Label label53;
            System.Windows.Forms.Label label54;
            System.Windows.Forms.Label label55;
            System.Windows.Forms.Label label56;
            System.Windows.Forms.Label label57;
            System.Windows.Forms.Label label58;
            System.Windows.Forms.Label label59;
            System.Windows.Forms.Label label64;
            System.Windows.Forms.Label label65;
            System.Windows.Forms.Label label66;
            System.Windows.Forms.Label label67;
            System.Windows.Forms.Label label60;
            System.Windows.Forms.Label label69;
            System.Windows.Forms.Label label70;
            System.Windows.Forms.Label eluviation_rate_constant;
            System.Windows.Forms.Label label72;
            System.Windows.Forms.Label label68;
            System.Windows.Forms.Label label71;
            System.Windows.Forms.Label label73;
            System.Windows.Forms.Label label74;
            System.Windows.Forms.Label label75;
            System.Windows.Forms.Label label76;
            System.Windows.Forms.Label label77;
            System.Windows.Forms.Label label81;
            System.Windows.Forms.Label label80;
            System.Windows.Forms.Label label82;
            System.Windows.Forms.Label label83;
            System.Windows.Forms.Label label84;
            System.Windows.Forms.Label label85;
            System.Windows.Forms.Label label86;
            System.Windows.Forms.Label label13;
            this.label36 = new System.Windows.Forms.Label();
            this.radio_ls_fraction = new System.Windows.Forms.RadioButton();
            this.radio_ls_absolute = new System.Windows.Forms.RadioButton();
            this.label35 = new System.Windows.Forms.Label();
            this.label34 = new System.Windows.Forms.Label();
            this.text_ls_rel_rain_intens = new System.Windows.Forms.TextBox();
            this.textBox_ls_trans = new System.Windows.Forms.TextBox();
            this.textBox_ls_bd = new System.Windows.Forms.TextBox();
            this.textBox_ls_ifr = new System.Windows.Forms.TextBox();
            this.textBox_ls_coh = new System.Windows.Forms.TextBox();
            this.text_ls_abs_rain_intens = new System.Windows.Forms.TextBox();
            this.label32 = new System.Windows.Forms.Label();
            this.label31 = new System.Windows.Forms.Label();
            this.label30 = new System.Windows.Forms.Label();
            this.label22 = new System.Windows.Forms.Label();
            this.label18 = new System.Windows.Forms.Label();
            this.pictureBox4 = new System.Windows.Forms.PictureBox();
            this.Landslide_checkbox = new System.Windows.Forms.CheckBox();
            this.mainMenu1 = new System.Windows.Forms.MainMenu(this.components);
            this.menuItemConfigFile = new System.Windows.Forms.MenuItem();
            this.menuItemConfigFileOpen = new System.Windows.Forms.MenuItem();
            this.menuItemConfigFileSaveAs = new System.Windows.Forms.MenuItem();
            this.menuItemConfigFileSave = new System.Windows.Forms.MenuItem();
            this.Mapselector = new System.Windows.Forms.MenuItem();
            this.Menu_map_total_sediment = new System.Windows.Forms.MenuItem();
            this.Menu_map_waterflow = new System.Windows.Forms.MenuItem();
            this.Menu_map_tillage = new System.Windows.Forms.MenuItem();
            this.Menu_map_water_ero = new System.Windows.Forms.MenuItem();
            this.Menu_map_creep = new System.Windows.Forms.MenuItem();
            this.Menu_map_weathering = new System.Windows.Forms.MenuItem();
            this.Menu_map_landsliding = new System.Windows.Forms.MenuItem();
            this.Menu_map_critical_rainfall = new System.Windows.Forms.MenuItem();
            this.Menu_About_box = new System.Windows.Forms.MenuItem();
            this.statusBar1 = new System.Windows.Forms.StatusBar();
            this.InfoStatusPanel = new System.Windows.Forms.StatusBarPanel();
            this.TimeStatusPanel = new System.Windows.Forms.StatusBarPanel();
            this.ProcessStatusPanel = new System.Windows.Forms.StatusBarPanel();
            this.out_sed_statuspanel = new System.Windows.Forms.StatusBarPanel();
            this.total_tillage_statuspanel = new System.Windows.Forms.StatusBarPanel();
            this.start_button = new System.Windows.Forms.Button();
            this.End_button = new System.Windows.Forms.Button();
            this.View_tabs_checkbox = new System.Windows.Forms.CheckBox();
            this.toolTip1 = new System.Windows.Forms.ToolTip(this.components);
            this.label2 = new System.Windows.Forms.Label();
            this.checkBoxGenerateAVIFile = new System.Windows.Forms.CheckBox();
            this.label33 = new System.Windows.Forms.Label();
            this.googleAnimationCheckbox = new System.Windows.Forms.CheckBox();
            this.label79 = new System.Windows.Forms.Label();
            this.textBoxAVIFile = new System.Windows.Forms.TextBox();
            this.googleAnimationTextBox = new System.Windows.Forms.TextBox();
            this.UTMzonebox = new System.Windows.Forms.TextBox();
            this.label25 = new System.Windows.Forms.Label();
            this.label14 = new System.Windows.Forms.Label();
            this.label15 = new System.Windows.Forms.Label();
            this.label17 = new System.Windows.Forms.Label();
            this.fill_sinks_before_checkbox = new System.Windows.Forms.CheckBox();
            this.label8 = new System.Windows.Forms.Label();
            this.fill_sinks_during_checkbox = new System.Windows.Forms.CheckBox();
            this.groupBox13 = new System.Windows.Forms.GroupBox();
            this.label7 = new System.Windows.Forms.Label();
            this.label5 = new System.Windows.Forms.Label();
            this.label29 = new System.Windows.Forms.Label();
            this.groupBox3 = new System.Windows.Forms.GroupBox();
            this.landuse_determinator_button = new System.Windows.Forms.Button();
            this.groupBox9 = new System.Windows.Forms.GroupBox();
            this.parameter_k1_textbox = new System.Windows.Forms.TextBox();
            this.label24 = new System.Windows.Forms.Label();
            this.label26 = new System.Windows.Forms.Label();
            this.label27 = new System.Windows.Forms.Label();
            this.label28 = new System.Windows.Forms.Label();
            this.parameter_k2_textbox = new System.Windows.Forms.TextBox();
            this.parameter_Pa_textbox = new System.Windows.Forms.TextBox();
            this.parameter_P0_textbox = new System.Windows.Forms.TextBox();
            this.label21 = new System.Windows.Forms.Label();
            this.Biological_weathering_checkbox = new System.Windows.Forms.CheckBox();
            this.label98 = new System.Windows.Forms.Label();
            this.label99 = new System.Windows.Forms.Label();
            this.trackBar1 = new System.Windows.Forms.TrackBar();
            this.label61 = new System.Windows.Forms.Label();
            this.map_controls = new System.Windows.Forms.GroupBox();
            this.label63 = new System.Windows.Forms.Label();
            this.graphicToGoogleEarthButton = new System.Windows.Forms.Button();
            this.trackBar2 = new System.Windows.Forms.TrackBar();
            this.label62 = new System.Windows.Forms.Label();
            this.comboBox1 = new System.Windows.Forms.ComboBox();
            this.openFileDialog1 = new System.Windows.Forms.OpenFileDialog();
            this.label1 = new System.Windows.Forms.Label();
            this.button6 = new System.Windows.Forms.Button();
            this.textBox1 = new System.Windows.Forms.TextBox();
            this.textBox2 = new System.Windows.Forms.TextBox();
            this.Output = new System.Windows.Forms.TabPage();
            this.groupBox6 = new System.Windows.Forms.GroupBox();
            this.groupBox12 = new System.Windows.Forms.GroupBox();
            this.annual_output_checkbox = new System.Windows.Forms.RadioButton();
            this.cumulative_output_checkbox = new System.Windows.Forms.RadioButton();
            this.groupBox11 = new System.Windows.Forms.GroupBox();
            this.Regular_output_checkbox = new System.Windows.Forms.CheckBox();
            this.Final_output_checkbox = new System.Windows.Forms.CheckBox();
            this.Box_years_output = new System.Windows.Forms.TextBox();
            this.groupBox1 = new System.Windows.Forms.GroupBox();
            this.diagnostic_output_checkbox = new System.Windows.Forms.CheckBox();
            this.label37 = new System.Windows.Forms.Label();
            this.outputcode_textbox = new System.Windows.Forms.TextBox();
            this.water_output_checkbox = new System.Windows.Forms.CheckBox();
            this.depressions_output_checkbox = new System.Windows.Forms.CheckBox();
            this.all_process_output_checkbox = new System.Windows.Forms.CheckBox();
            this.Soildepth_output_checkbox = new System.Windows.Forms.CheckBox();
            this.Alt_change_output_checkbox = new System.Windows.Forms.CheckBox();
            this.Altitude_output_checkbox = new System.Windows.Forms.CheckBox();
            this.checkedListBox1 = new System.Windows.Forms.CheckedListBox();
            this.groupBox5 = new System.Windows.Forms.GroupBox();
            this.button3 = new System.Windows.Forms.Button();
            this.button1 = new System.Windows.Forms.Button();
            this.button2 = new System.Windows.Forms.Button();
            this.timeseries_form_button = new System.Windows.Forms.Button();
            this.UTMgroupBox = new System.Windows.Forms.GroupBox();
            this.textBox6 = new System.Windows.Forms.TextBox();
            this.UTMsouthcheck = new System.Windows.Forms.CheckBox();
            this.UTMgridcheckbox = new System.Windows.Forms.CheckBox();
            this.textBox4 = new System.Windows.Forms.TextBox();
            this.googleBeginDate = new System.Windows.Forms.TextBox();
            this.googAnimationSaveInterval = new System.Windows.Forms.TextBox();
            this.saveintervalbox = new System.Windows.Forms.TextBox();
            this.label78 = new System.Windows.Forms.Label();
            this.Run = new System.Windows.Forms.TabPage();
            this.button4 = new System.Windows.Forms.Button();
            this.version_lux_checkbox = new System.Windows.Forms.CheckBox();
            this.groupBox2 = new System.Windows.Forms.GroupBox();
            this.label120 = new System.Windows.Forms.Label();
            this.calibration_ratio_reduction_parameter_textbox = new System.Windows.Forms.TextBox();
            this.label119 = new System.Windows.Forms.Label();
            this.calibration_levels_textbox = new System.Windows.Forms.TextBox();
            this.label116 = new System.Windows.Forms.Label();
            this.label118 = new System.Windows.Forms.Label();
            this.label117 = new System.Windows.Forms.Label();
            this.label115 = new System.Windows.Forms.Label();
            this.label114 = new System.Windows.Forms.Label();
            this.Sensitivity_button = new System.Windows.Forms.RadioButton();
            this.Calibration_button = new System.Windows.Forms.RadioButton();
            this.label113 = new System.Windows.Forms.Label();
            this.calibration_ratios_textbox = new System.Windows.Forms.TextBox();
            this.calibration = new System.Windows.Forms.CheckBox();
            this.Ik_ben_Marijn = new System.Windows.Forms.CheckBox();
            this.groupBox7 = new System.Windows.Forms.GroupBox();
            this.runs_checkbox = new System.Windows.Forms.RadioButton();
            this.label16 = new System.Windows.Forms.Label();
            this.Number_runs_textbox = new System.Windows.Forms.TextBox();
            this.Input = new System.Windows.Forms.TabPage();
            this.check_time_T = new System.Windows.Forms.CheckBox();
            this.temp_input_filename_textbox = new System.Windows.Forms.TextBox();
            this.temp_constant_value_box = new System.Windows.Forms.TextBox();
            this.soil_specify_button = new System.Windows.Forms.Button();
            this.label88 = new System.Windows.Forms.Label();
            this.explain_input_button = new System.Windows.Forms.Button();
            this.check_time_evap = new System.Windows.Forms.CheckBox();
            this.check_time_infil = new System.Windows.Forms.CheckBox();
            this.check_time_rain = new System.Windows.Forms.CheckBox();
            this.check_time_till_fields = new System.Windows.Forms.CheckBox();
            this.check_time_landuse = new System.Windows.Forms.CheckBox();
            this.check_space_DTM = new System.Windows.Forms.CheckBox();
            this.tillfields_constant_textbox = new System.Windows.Forms.TextBox();
            this.tillfields_input_filename_textbox = new System.Windows.Forms.TextBox();
            this.evap_constant_value_box = new System.Windows.Forms.TextBox();
            this.evap_input_filename_textbox = new System.Windows.Forms.TextBox();
            this.infil_constant_value_box = new System.Windows.Forms.TextBox();
            this.infil_input_filename_textbox = new System.Windows.Forms.TextBox();
            this.rainfall_constant_value_box = new System.Windows.Forms.TextBox();
            this.landuse_constant_value_box = new System.Windows.Forms.TextBox();
            this.soildepth_constant_value_box = new System.Windows.Forms.TextBox();
            this.landuse_input_filename_textbox = new System.Windows.Forms.TextBox();
            this.soildepth_input_filename_textbox = new System.Windows.Forms.TextBox();
            this.rain_input_filename_textbox = new System.Windows.Forms.TextBox();
            this.dtm_input_filename_textbox = new System.Windows.Forms.TextBox();
            this.groupBox8 = new System.Windows.Forms.GroupBox();
            this.check_space_evap = new System.Windows.Forms.CheckBox();
            this.check_space_infil = new System.Windows.Forms.CheckBox();
            this.check_space_rain = new System.Windows.Forms.CheckBox();
            this.check_space_till_fields = new System.Windows.Forms.CheckBox();
            this.check_space_landuse = new System.Windows.Forms.CheckBox();
            this.check_space_soildepth = new System.Windows.Forms.CheckBox();
            this.label4 = new System.Windows.Forms.Label();
            this.label3 = new System.Windows.Forms.Label();
            this.label23 = new System.Windows.Forms.Label();
            this.Processes = new System.Windows.Forms.TabPage();
            this.Process_tabs = new System.Windows.Forms.TabControl();
            this.Water = new System.Windows.Forms.TabPage();
            this.daily_water = new System.Windows.Forms.CheckBox();
            this.label87 = new System.Windows.Forms.Label();
            this.selectivity_constant_textbox = new System.Windows.Forms.TextBox();
            this.bio_protection_constant_textbox = new System.Windows.Forms.TextBox();
            this.erosion_threshold_textbox = new System.Windows.Forms.TextBox();
            this.rock_protection_constant_textbox = new System.Windows.Forms.TextBox();
            this.label90 = new System.Windows.Forms.Label();
            this.label91 = new System.Windows.Forms.Label();
            this.label92 = new System.Windows.Forms.Label();
            this.parameter_n_textbox = new System.Windows.Forms.TextBox();
            this.parameter_conv_textbox = new System.Windows.Forms.TextBox();
            this.parameter_K_textbox = new System.Windows.Forms.TextBox();
            this.parameter_m_textbox = new System.Windows.Forms.TextBox();
            this.only_waterflow_checkbox = new System.Windows.Forms.CheckBox();
            this.pictureBox1 = new System.Windows.Forms.PictureBox();
            this.label12 = new System.Windows.Forms.Label();
            this.label11 = new System.Windows.Forms.Label();
            this.label10 = new System.Windows.Forms.Label();
            this.label9 = new System.Windows.Forms.Label();
            this.Water_ero_checkbox = new System.Windows.Forms.CheckBox();
            this.Tillage = new System.Windows.Forms.TabPage();
            this.pictureBox2 = new System.Windows.Forms.PictureBox();
            this.label20 = new System.Windows.Forms.Label();
            this.trte = new System.Windows.Forms.Label();
            this.parameter_tillage_constant_textbox = new System.Windows.Forms.TextBox();
            this.parameter_ploughing_depth_textbox = new System.Windows.Forms.TextBox();
            this.Tillage_checkbox = new System.Windows.Forms.CheckBox();
            this.Creeper = new System.Windows.Forms.TabPage();
            this.creep_testing = new System.Windows.Forms.CheckBox();
            this.pictureBox3 = new System.Windows.Forms.PictureBox();
            this.label19 = new System.Windows.Forms.Label();
            this.parameter_diffusivity_textbox = new System.Windows.Forms.TextBox();
            this.creep_active_checkbox = new System.Windows.Forms.CheckBox();
            this.Solifluction = new System.Windows.Forms.TabPage();
            this.pictureBox5 = new System.Windows.Forms.PictureBox();
            this.Solifluction_checkbox = new System.Windows.Forms.CheckBox();
            this.Rock_weathering = new System.Windows.Forms.TabPage();
            this.rockweath_method = new System.Windows.Forms.ComboBox();
            this.pictureBox6 = new System.Windows.Forms.PictureBox();
            this.groupBox10 = new System.Windows.Forms.GroupBox();
            this.Frost_weathering_checkbox = new System.Windows.Forms.CheckBox();
            this.Tectonics = new System.Windows.Forms.TabPage();
            this.groupBox14 = new System.Windows.Forms.GroupBox();
            this.groupBox16 = new System.Windows.Forms.GroupBox();
            this.text_lift_col_less = new System.Windows.Forms.TextBox();
            this.text_lift_col_more = new System.Windows.Forms.TextBox();
            this.text_lift_row_less = new System.Windows.Forms.TextBox();
            this.text_lift_row_more = new System.Windows.Forms.TextBox();
            this.radio_lift_col_less_than = new System.Windows.Forms.RadioButton();
            this.radio_lift_row_more_than = new System.Windows.Forms.RadioButton();
            this.radio_lift_col_more_than = new System.Windows.Forms.RadioButton();
            this.radio_lift_row_less_than = new System.Windows.Forms.RadioButton();
            this.Uplift_rate_textbox = new System.Windows.Forms.TextBox();
            this.uplift_active_checkbox = new System.Windows.Forms.CheckBox();
            this.label39 = new System.Windows.Forms.Label();
            this.groupBox4 = new System.Windows.Forms.GroupBox();
            this.label38 = new System.Windows.Forms.Label();
            this.Tilting_rate_textbox = new System.Windows.Forms.TextBox();
            this.groupBox15 = new System.Windows.Forms.GroupBox();
            this.radio_tilt_col_max = new System.Windows.Forms.RadioButton();
            this.radio_tilt_row_zero = new System.Windows.Forms.RadioButton();
            this.radio_tilt_col_zero = new System.Windows.Forms.RadioButton();
            this.radio_tilt_row_max = new System.Windows.Forms.RadioButton();
            this.tilting_active_checkbox = new System.Windows.Forms.CheckBox();
            this.treefall = new System.Windows.Forms.TabPage();
            this.tf_freq = new System.Windows.Forms.TextBox();
            this.label112 = new System.Windows.Forms.Label();
            this.tf_age = new System.Windows.Forms.TextBox();
            this.label111 = new System.Windows.Forms.Label();
            this.tf_growth = new System.Windows.Forms.TextBox();
            this.label110 = new System.Windows.Forms.Label();
            this.tf_D = new System.Windows.Forms.TextBox();
            this.label95 = new System.Windows.Forms.Label();
            this.label107 = new System.Windows.Forms.Label();
            this.tf_W = new System.Windows.Forms.TextBox();
            this.treefall_checkbox = new System.Windows.Forms.CheckBox();
            this.Creep_Checkbox = new System.Windows.Forms.CheckBox();
            this.tabControl1 = new System.Windows.Forms.TabControl();
            this.tabPage1 = new System.Windows.Forms.TabPage();
            this.tabControl2 = new System.Windows.Forms.TabControl();
            this.physical = new System.Windows.Forms.TabPage();
            this.upper_particle_fine_clay_textbox = new System.Windows.Forms.TextBox();
            this.upper_particle_clay_textbox = new System.Windows.Forms.TextBox();
            this.upper_particle_silt_textbox = new System.Windows.Forms.TextBox();
            this.upper_particle_sand_textbox = new System.Windows.Forms.TextBox();
            this.upper_particle_coarse_textbox = new System.Windows.Forms.TextBox();
            this.physical_weath_constant2 = new System.Windows.Forms.TextBox();
            this.physical_weath_constant1 = new System.Windows.Forms.TextBox();
            this.Physical_weath_C1_textbox = new System.Windows.Forms.TextBox();
            this.soil_phys_weath_checkbox = new System.Windows.Forms.CheckBox();
            this.chemical = new System.Windows.Forms.TabPage();
            this.specific_area_fine_clay_textbox = new System.Windows.Forms.TextBox();
            this.specific_area_clay_textbox = new System.Windows.Forms.TextBox();
            this.specific_area_silt_textbox = new System.Windows.Forms.TextBox();
            this.specific_area_sand_textbox = new System.Windows.Forms.TextBox();
            this.specific_area_coarse_textbox = new System.Windows.Forms.TextBox();
            this.chem_weath_specific_coefficient_textbox = new System.Windows.Forms.TextBox();
            this.chem_weath_depth_constant_textbox = new System.Windows.Forms.TextBox();
            this.chem_weath_rate_constant_textbox = new System.Windows.Forms.TextBox();
            this.soil_chem_weath_checkbox = new System.Windows.Forms.CheckBox();
            this.clay = new System.Windows.Forms.TabPage();
            this.ct_Jagercikova = new System.Windows.Forms.CheckBox();
            this.label109 = new System.Windows.Forms.Label();
            this.label108 = new System.Windows.Forms.Label();
            this.ct_dd_Jagercikova = new System.Windows.Forms.TextBox();
            this.ct_v0_Jagercikova = new System.Windows.Forms.TextBox();
            this.ct_depth_decay = new System.Windows.Forms.TextBox();
            this.CT_depth_decay_checkbox = new System.Windows.Forms.CheckBox();
            this.eluviation_coefficient_textbox = new System.Windows.Forms.TextBox();
            this.maximum_eluviation_textbox = new System.Windows.Forms.TextBox();
            this.clay_neoform_C2_textbox = new System.Windows.Forms.TextBox();
            this.clay_neoform_C1_textbox = new System.Windows.Forms.TextBox();
            this.clay_neoform_constant_textbox = new System.Windows.Forms.TextBox();
            this.soil_clay_transloc_checkbox = new System.Windows.Forms.CheckBox();
            this.bioturbation = new System.Windows.Forms.TabPage();
            this.bioturbation_depth_decay_textbox = new System.Windows.Forms.TextBox();
            this.potential_bioturbation_textbox = new System.Windows.Forms.TextBox();
            this.soil_bioturb_checkbox = new System.Windows.Forms.CheckBox();
            this.carbon = new System.Windows.Forms.TabPage();
            this.carbon_o_decomp_rate_textbox = new System.Windows.Forms.TextBox();
            this.carbon_y_decomp_rate_textbox = new System.Windows.Forms.TextBox();
            this.carbon_o_twi_decay_textbox = new System.Windows.Forms.TextBox();
            this.carbon_y_twi_decay_textbox = new System.Windows.Forms.TextBox();
            this.carbon_o_depth_decay_textbox = new System.Windows.Forms.TextBox();
            this.carbon_y_depth_decay_textbox = new System.Windows.Forms.TextBox();
            this.carbon_humification_fraction_textbox = new System.Windows.Forms.TextBox();
            this.carbon_depth_decay_textbox = new System.Windows.Forms.TextBox();
            this.carbon_input_textbox = new System.Windows.Forms.TextBox();
            this.soil_carbon_cycle_checkbox = new System.Windows.Forms.CheckBox();
            this.decalcification = new System.Windows.Forms.TabPage();
            this.label94 = new System.Windows.Forms.Label();
            this.ini_CaCO3_content = new System.Windows.Forms.TextBox();
            this.decalcification_checkbox = new System.Windows.Forms.CheckBox();
            this.tabPage2 = new System.Windows.Forms.TabPage();
            this.check_scaling_daily_weather = new System.Windows.Forms.CheckBox();
            this.label106 = new System.Windows.Forms.Label();
            this.snow_threshold_textbox = new System.Windows.Forms.TextBox();
            this.label105 = new System.Windows.Forms.Label();
            this.snowmelt_factor_textbox = new System.Windows.Forms.TextBox();
            this.label104 = new System.Windows.Forms.Label();
            this.latitude_min = new System.Windows.Forms.TextBox();
            this.label103 = new System.Windows.Forms.Label();
            this.latitude_deg = new System.Windows.Forms.TextBox();
            this.label100 = new System.Windows.Forms.Label();
            this.label101 = new System.Windows.Forms.Label();
            this.label102 = new System.Windows.Forms.Label();
            this.dailyT_min = new System.Windows.Forms.TextBox();
            this.dailyT_max = new System.Windows.Forms.TextBox();
            this.dailyT_avg = new System.Windows.Forms.TextBox();
            this.label97 = new System.Windows.Forms.Label();
            this.daily_n = new System.Windows.Forms.TextBox();
            this.label96 = new System.Windows.Forms.Label();
            this.label93 = new System.Windows.Forms.Label();
            this.label89 = new System.Windows.Forms.Label();
            this.label40 = new System.Windows.Forms.Label();
            this.dailyET0 = new System.Windows.Forms.TextBox();
            this.dailyD = new System.Windows.Forms.TextBox();
            this.dailyP = new System.Windows.Forms.TextBox();
            this.view_maps_checkbox = new System.Windows.Forms.CheckBox();
            this.timer1 = new System.Windows.Forms.Timer(this.components);
            label6 = new System.Windows.Forms.Label();
            Landsliding = new System.Windows.Forms.TabPage();
            label41 = new System.Windows.Forms.Label();
            label42 = new System.Windows.Forms.Label();
            label43 = new System.Windows.Forms.Label();
            label44 = new System.Windows.Forms.Label();
            label45 = new System.Windows.Forms.Label();
            label46 = new System.Windows.Forms.Label();
            label47 = new System.Windows.Forms.Label();
            label48 = new System.Windows.Forms.Label();
            label49 = new System.Windows.Forms.Label();
            label50 = new System.Windows.Forms.Label();
            label51 = new System.Windows.Forms.Label();
            label52 = new System.Windows.Forms.Label();
            label53 = new System.Windows.Forms.Label();
            label54 = new System.Windows.Forms.Label();
            label55 = new System.Windows.Forms.Label();
            label56 = new System.Windows.Forms.Label();
            label57 = new System.Windows.Forms.Label();
            label58 = new System.Windows.Forms.Label();
            label59 = new System.Windows.Forms.Label();
            label64 = new System.Windows.Forms.Label();
            label65 = new System.Windows.Forms.Label();
            label66 = new System.Windows.Forms.Label();
            label67 = new System.Windows.Forms.Label();
            label60 = new System.Windows.Forms.Label();
            label69 = new System.Windows.Forms.Label();
            label70 = new System.Windows.Forms.Label();
            eluviation_rate_constant = new System.Windows.Forms.Label();
            label72 = new System.Windows.Forms.Label();
            label68 = new System.Windows.Forms.Label();
            label71 = new System.Windows.Forms.Label();
            label73 = new System.Windows.Forms.Label();
            label74 = new System.Windows.Forms.Label();
            label75 = new System.Windows.Forms.Label();
            label76 = new System.Windows.Forms.Label();
            label77 = new System.Windows.Forms.Label();
            label81 = new System.Windows.Forms.Label();
            label80 = new System.Windows.Forms.Label();
            label82 = new System.Windows.Forms.Label();
            label83 = new System.Windows.Forms.Label();
            label84 = new System.Windows.Forms.Label();
            label85 = new System.Windows.Forms.Label();
            label86 = new System.Windows.Forms.Label();
            label13 = new System.Windows.Forms.Label();
            Landsliding.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox4)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.InfoStatusPanel)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.TimeStatusPanel)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.ProcessStatusPanel)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.out_sed_statuspanel)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.total_tillage_statuspanel)).BeginInit();
            this.groupBox13.SuspendLayout();
            this.groupBox3.SuspendLayout();
            this.groupBox9.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.trackBar1)).BeginInit();
            this.map_controls.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.trackBar2)).BeginInit();
            this.Output.SuspendLayout();
            this.groupBox6.SuspendLayout();
            this.groupBox12.SuspendLayout();
            this.groupBox11.SuspendLayout();
            this.groupBox1.SuspendLayout();
            this.groupBox5.SuspendLayout();
            this.UTMgroupBox.SuspendLayout();
            this.Run.SuspendLayout();
            this.groupBox2.SuspendLayout();
            this.groupBox7.SuspendLayout();
            this.Input.SuspendLayout();
            this.groupBox8.SuspendLayout();
            this.Processes.SuspendLayout();
            this.Process_tabs.SuspendLayout();
            this.Water.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox1)).BeginInit();
            this.Tillage.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox2)).BeginInit();
            this.Creeper.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox3)).BeginInit();
            this.Solifluction.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox5)).BeginInit();
            this.Rock_weathering.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox6)).BeginInit();
            this.groupBox10.SuspendLayout();
            this.Tectonics.SuspendLayout();
            this.groupBox14.SuspendLayout();
            this.groupBox16.SuspendLayout();
            this.groupBox4.SuspendLayout();
            this.groupBox15.SuspendLayout();
            this.treefall.SuspendLayout();
            this.tabControl1.SuspendLayout();
            this.tabPage1.SuspendLayout();
            this.tabControl2.SuspendLayout();
            this.physical.SuspendLayout();
            this.chemical.SuspendLayout();
            this.clay.SuspendLayout();
            this.bioturbation.SuspendLayout();
            this.carbon.SuspendLayout();
            this.decalcification.SuspendLayout();
            this.tabPage2.SuspendLayout();
            this.SuspendLayout();
            // 
            // label6
            // 
            label6.Location = new System.Drawing.Point(179, 22);
            label6.Name = "label6";
            label6.Size = new System.Drawing.Size(41, 24);
            label6.TabIndex = 109;
            label6.Text = "f(x,y)";
            label6.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            this.toolTip1.SetToolTip(label6, "Check this for spatially variable inputs");
            // 
            // Landsliding
            // 
            Landsliding.Controls.Add(this.label36);
            Landsliding.Controls.Add(this.radio_ls_fraction);
            Landsliding.Controls.Add(this.radio_ls_absolute);
            Landsliding.Controls.Add(this.label35);
            Landsliding.Controls.Add(this.label34);
            Landsliding.Controls.Add(this.text_ls_rel_rain_intens);
            Landsliding.Controls.Add(this.textBox_ls_trans);
            Landsliding.Controls.Add(this.textBox_ls_bd);
            Landsliding.Controls.Add(this.textBox_ls_ifr);
            Landsliding.Controls.Add(this.textBox_ls_coh);
            Landsliding.Controls.Add(this.text_ls_abs_rain_intens);
            Landsliding.Controls.Add(this.label32);
            Landsliding.Controls.Add(this.label31);
            Landsliding.Controls.Add(this.label30);
            Landsliding.Controls.Add(this.label22);
            Landsliding.Controls.Add(this.label18);
            Landsliding.Controls.Add(this.pictureBox4);
            Landsliding.Controls.Add(this.Landslide_checkbox);
            Landsliding.Location = new System.Drawing.Point(4, 22);
            Landsliding.Name = "Landsliding";
            Landsliding.Size = new System.Drawing.Size(732, 250);
            Landsliding.TabIndex = 2;
            Landsliding.Text = "Landsliding";
            Landsliding.UseVisualStyleBackColor = true;
            // 
            // label36
            // 
            this.label36.AutoSize = true;
            this.label36.Location = new System.Drawing.Point(49, 115);
            this.label36.Name = "label36";
            this.label36.Size = new System.Drawing.Size(236, 13);
            this.label36.TabIndex = 30;
            this.label36.Text = "Parameters for critical rainfall intensity calculation";
            // 
            // radio_ls_fraction
            // 
            this.radio_ls_fraction.AutoSize = true;
            this.radio_ls_fraction.Checked = true;
            this.radio_ls_fraction.Location = new System.Drawing.Point(53, 83);
            this.radio_ls_fraction.Name = "radio_ls_fraction";
            this.radio_ls_fraction.Size = new System.Drawing.Size(14, 13);
            this.radio_ls_fraction.TabIndex = 29;
            this.radio_ls_fraction.TabStop = true;
            this.radio_ls_fraction.UseVisualStyleBackColor = true;
            // 
            // radio_ls_absolute
            // 
            this.radio_ls_absolute.AutoSize = true;
            this.radio_ls_absolute.Location = new System.Drawing.Point(53, 57);
            this.radio_ls_absolute.Name = "radio_ls_absolute";
            this.radio_ls_absolute.Size = new System.Drawing.Size(14, 13);
            this.radio_ls_absolute.TabIndex = 28;
            this.radio_ls_absolute.UseVisualStyleBackColor = true;
            // 
            // label35
            // 
            this.label35.AutoSize = true;
            this.label35.Enabled = false;
            this.label35.Location = new System.Drawing.Point(138, 57);
            this.label35.Name = "label35";
            this.label35.Size = new System.Drawing.Size(105, 13);
            this.label35.TabIndex = 27;
            this.label35.Text = "Absolute value [m/d]";
            // 
            // label34
            // 
            this.label34.AutoSize = true;
            this.label34.Location = new System.Drawing.Point(138, 83);
            this.label34.Name = "label34";
            this.label34.Size = new System.Drawing.Size(270, 13);
            this.label34.TabIndex = 26;
            this.label34.Text = "Fraction of total annual rainfall [between 1 and 0.00274]";
            // 
            // text_ls_rel_rain_intens
            // 
            this.text_ls_rel_rain_intens.Location = new System.Drawing.Point(79, 80);
            this.text_ls_rel_rain_intens.Name = "text_ls_rel_rain_intens";
            this.text_ls_rel_rain_intens.Size = new System.Drawing.Size(53, 20);
            this.text_ls_rel_rain_intens.TabIndex = 25;
            this.text_ls_rel_rain_intens.Text = "0.1";
            // 
            // textBox_ls_trans
            // 
            this.textBox_ls_trans.Location = new System.Drawing.Point(52, 210);
            this.textBox_ls_trans.Name = "textBox_ls_trans";
            this.textBox_ls_trans.Size = new System.Drawing.Size(53, 20);
            this.textBox_ls_trans.TabIndex = 23;
            this.textBox_ls_trans.Text = "15";
            // 
            // textBox_ls_bd
            // 
            this.textBox_ls_bd.Location = new System.Drawing.Point(52, 184);
            this.textBox_ls_bd.Name = "textBox_ls_bd";
            this.textBox_ls_bd.Size = new System.Drawing.Size(53, 20);
            this.textBox_ls_bd.TabIndex = 21;
            this.textBox_ls_bd.Text = "1.4";
            // 
            // textBox_ls_ifr
            // 
            this.textBox_ls_ifr.Location = new System.Drawing.Point(52, 158);
            this.textBox_ls_ifr.Name = "textBox_ls_ifr";
            this.textBox_ls_ifr.Size = new System.Drawing.Size(53, 20);
            this.textBox_ls_ifr.TabIndex = 19;
            this.textBox_ls_ifr.Text = "0.7";
            // 
            // textBox_ls_coh
            // 
            this.textBox_ls_coh.Location = new System.Drawing.Point(52, 131);
            this.textBox_ls_coh.Name = "textBox_ls_coh";
            this.textBox_ls_coh.Size = new System.Drawing.Size(53, 20);
            this.textBox_ls_coh.TabIndex = 17;
            this.textBox_ls_coh.Text = "0.15";
            // 
            // text_ls_abs_rain_intens
            // 
            this.text_ls_abs_rain_intens.Enabled = false;
            this.text_ls_abs_rain_intens.Location = new System.Drawing.Point(79, 54);
            this.text_ls_abs_rain_intens.Name = "text_ls_abs_rain_intens";
            this.text_ls_abs_rain_intens.Size = new System.Drawing.Size(53, 20);
            this.text_ls_abs_rain_intens.TabIndex = 15;
            this.text_ls_abs_rain_intens.Text = "0.1";
            // 
            // label32
            // 
            this.label32.AutoSize = true;
            this.label32.Location = new System.Drawing.Point(111, 213);
            this.label32.Name = "label32";
            this.label32.Size = new System.Drawing.Size(169, 13);
            this.label32.TabIndex = 24;
            this.label32.Text = "Saturated soil transmissivity [m2/d]";
            // 
            // label31
            // 
            this.label31.AutoSize = true;
            this.label31.Location = new System.Drawing.Point(111, 187);
            this.label31.Name = "label31";
            this.label31.Size = new System.Drawing.Size(105, 13);
            this.label31.TabIndex = 22;
            this.label31.Text = "Bulk density [kg m-3]";
            // 
            // label30
            // 
            this.label30.AutoSize = true;
            this.label30.Location = new System.Drawing.Point(111, 161);
            this.label30.Name = "label30";
            this.label30.Size = new System.Drawing.Size(152, 13);
            this.label30.TabIndex = 20;
            this.label30.Text = "Internal friction angle [degrees]";
            // 
            // label22
            // 
            this.label22.AutoSize = true;
            this.label22.Location = new System.Drawing.Point(111, 134);
            this.label22.Name = "label22";
            this.label22.Size = new System.Drawing.Size(112, 13);
            this.label22.TabIndex = 18;
            this.label22.Text = "Combined cohesion [-]";
            // 
            // label18
            // 
            this.label18.AutoSize = true;
            this.label18.Location = new System.Drawing.Point(50, 38);
            this.label18.Name = "label18";
            this.label18.Size = new System.Drawing.Size(117, 13);
            this.label18.TabIndex = 16;
            this.label18.Text = "Critical rainfall threshold";
            // 
            // pictureBox4
            // 
            this.pictureBox4.Image = ((System.Drawing.Image)(resources.GetObject("pictureBox4.Image")));
            this.pictureBox4.Location = new System.Drawing.Point(480, 57);
            this.pictureBox4.Name = "pictureBox4";
            this.pictureBox4.Size = new System.Drawing.Size(180, 137);
            this.pictureBox4.TabIndex = 14;
            this.pictureBox4.TabStop = false;
            // 
            // Landslide_checkbox
            // 
            this.Landslide_checkbox.AutoSize = true;
            this.Landslide_checkbox.Location = new System.Drawing.Point(26, 14);
            this.Landslide_checkbox.Name = "Landslide_checkbox";
            this.Landslide_checkbox.Size = new System.Drawing.Size(124, 17);
            this.Landslide_checkbox.TabIndex = 1;
            this.Landslide_checkbox.Text = "Activate this process";
            this.Landslide_checkbox.UseVisualStyleBackColor = true;
            // 
            // label41
            // 
            label41.AutoSize = true;
            label41.Location = new System.Drawing.Point(142, 49);
            label41.Name = "label41";
            label41.Size = new System.Drawing.Size(147, 13);
            label41.TabIndex = 10;
            label41.Text = "weathering rate constant [y-1]";
            // 
            // label42
            // 
            label42.AutoSize = true;
            label42.Location = new System.Drawing.Point(142, 72);
            label42.Name = "label42";
            label42.Size = new System.Drawing.Size(136, 13);
            label42.TabIndex = 11;
            label42.Text = "depth decay constant [m-1]";
            // 
            // label43
            // 
            label43.AutoSize = true;
            label43.Location = new System.Drawing.Point(142, 98);
            label43.Name = "label43";
            label43.Size = new System.Drawing.Size(123, 13);
            label43.TabIndex = 12;
            label43.Text = "particle size constant [m]";
            // 
            // label44
            // 
            label44.AutoSize = true;
            label44.Location = new System.Drawing.Point(409, 49);
            label44.Name = "label44";
            label44.Size = new System.Drawing.Size(77, 13);
            label44.TabIndex = 13;
            label44.Text = "coarse fraction";
            // 
            // label45
            // 
            label45.AutoSize = true;
            label45.Location = new System.Drawing.Point(409, 72);
            label45.Name = "label45";
            label45.Size = new System.Drawing.Size(68, 13);
            label45.TabIndex = 14;
            label45.Text = "sand fraction";
            // 
            // label46
            // 
            label46.AutoSize = true;
            label46.Location = new System.Drawing.Point(409, 98);
            label46.Name = "label46";
            label46.Size = new System.Drawing.Size(57, 13);
            label46.TabIndex = 15;
            label46.Text = "silt fraction";
            // 
            // label47
            // 
            label47.AutoSize = true;
            label47.Location = new System.Drawing.Point(409, 124);
            label47.Name = "label47";
            label47.Size = new System.Drawing.Size(64, 13);
            label47.TabIndex = 16;
            label47.Text = "clay fraction";
            // 
            // label48
            // 
            label48.AutoSize = true;
            label48.Location = new System.Drawing.Point(409, 150);
            label48.Name = "label48";
            label48.Size = new System.Drawing.Size(84, 13);
            label48.TabIndex = 17;
            label48.Text = "fine clay fraction";
            // 
            // label49
            // 
            label49.AutoSize = true;
            label49.Location = new System.Drawing.Point(300, 27);
            label49.Name = "label49";
            label49.Size = new System.Drawing.Size(229, 13);
            label49.TabIndex = 18;
            label49.Text = "upper limit of particle size for texture classes [m]";
            // 
            // label50
            // 
            label50.AutoSize = true;
            label50.Location = new System.Drawing.Point(136, 93);
            label50.Name = "label50";
            label50.Size = new System.Drawing.Size(0, 13);
            label50.TabIndex = 18;
            // 
            // label51
            // 
            label51.AutoSize = true;
            label51.Location = new System.Drawing.Point(136, 67);
            label51.Name = "label51";
            label51.Size = new System.Drawing.Size(136, 13);
            label51.TabIndex = 17;
            label51.Text = "depth decay constant [m-1]";
            // 
            // label52
            // 
            label52.AutoSize = true;
            label52.Location = new System.Drawing.Point(136, 38);
            label52.Name = "label52";
            label52.Size = new System.Drawing.Size(268, 13);
            label52.TabIndex = 16;
            label52.Text = "weathering rate constant [kg / m2 mineral surface area]";
            // 
            // label53
            // 
            label53.AutoSize = true;
            label53.Location = new System.Drawing.Point(135, 97);
            label53.Name = "label53";
            label53.Size = new System.Drawing.Size(131, 13);
            label53.TabIndex = 19;
            label53.Text = "specific area coefficient [-]";
            // 
            // label54
            // 
            label54.AutoSize = true;
            label54.Location = new System.Drawing.Point(417, 22);
            label54.Name = "label54";
            label54.Size = new System.Drawing.Size(239, 13);
            label54.TabIndex = 30;
            label54.Text = "specific surface area for texture classes [m2 / kg]";
            // 
            // label55
            // 
            label55.AutoSize = true;
            label55.Location = new System.Drawing.Point(526, 145);
            label55.Name = "label55";
            label55.Size = new System.Drawing.Size(84, 13);
            label55.TabIndex = 29;
            label55.Text = "fine clay fraction";
            // 
            // label56
            // 
            label56.AutoSize = true;
            label56.Location = new System.Drawing.Point(526, 119);
            label56.Name = "label56";
            label56.Size = new System.Drawing.Size(64, 13);
            label56.TabIndex = 28;
            label56.Text = "clay fraction";
            // 
            // label57
            // 
            label57.AutoSize = true;
            label57.Location = new System.Drawing.Point(526, 93);
            label57.Name = "label57";
            label57.Size = new System.Drawing.Size(57, 13);
            label57.TabIndex = 27;
            label57.Text = "silt fraction";
            // 
            // label58
            // 
            label58.AutoSize = true;
            label58.Location = new System.Drawing.Point(526, 67);
            label58.Name = "label58";
            label58.Size = new System.Drawing.Size(68, 13);
            label58.TabIndex = 26;
            label58.Text = "sand fraction";
            // 
            // label59
            // 
            label59.AutoSize = true;
            label59.Location = new System.Drawing.Point(526, 44);
            label59.Name = "label59";
            label59.Size = new System.Drawing.Size(77, 13);
            label59.TabIndex = 25;
            label59.Text = "coarse fraction";
            // 
            // label64
            // 
            label64.AutoSize = true;
            label64.Location = new System.Drawing.Point(131, 134);
            label64.Name = "label64";
            label64.Size = new System.Drawing.Size(83, 13);
            label64.TabIndex = 46;
            label64.Text = "constant 2 [m-1]";
            // 
            // label65
            // 
            label65.AutoSize = true;
            label65.Location = new System.Drawing.Point(132, 130);
            label65.Name = "label65";
            label65.Size = new System.Drawing.Size(0, 13);
            label65.TabIndex = 45;
            // 
            // label66
            // 
            label66.AutoSize = true;
            label66.Location = new System.Drawing.Point(132, 104);
            label66.Name = "label66";
            label66.Size = new System.Drawing.Size(60, 13);
            label66.TabIndex = 44;
            label66.Text = "constant 1 ";
            // 
            // label67
            // 
            label67.AutoSize = true;
            label67.Location = new System.Drawing.Point(132, 75);
            label67.Name = "label67";
            label67.Size = new System.Drawing.Size(121, 13);
            label67.TabIndex = 43;
            label67.Text = "neoformation constant []";
            // 
            // label60
            // 
            label60.AutoSize = true;
            label60.Location = new System.Drawing.Point(23, 59);
            label60.Name = "label60";
            label60.Size = new System.Drawing.Size(110, 13);
            label60.TabIndex = 39;
            label60.Text = "fine clay neoformation";
            // 
            // label69
            // 
            label69.AutoSize = true;
            label69.Location = new System.Drawing.Point(411, 130);
            label69.Name = "label69";
            label69.Size = new System.Drawing.Size(0, 13);
            label69.TabIndex = 53;
            // 
            // label70
            // 
            label70.AutoSize = true;
            label70.Location = new System.Drawing.Point(411, 104);
            label70.Name = "label70";
            label70.Size = new System.Drawing.Size(97, 13);
            label70.TabIndex = 52;
            label70.Text = "saturation constant";
            // 
            // eluviation_rate_constant
            // 
            eluviation_rate_constant.AutoSize = true;
            eluviation_rate_constant.Location = new System.Drawing.Point(411, 75);
            eluviation_rate_constant.Name = "eluviation_rate_constant";
            eluviation_rate_constant.Size = new System.Drawing.Size(119, 13);
            eluviation_rate_constant.TabIndex = 51;
            eluviation_rate_constant.Text = "maximum eluviation [kg]";
            // 
            // label72
            // 
            label72.AutoSize = true;
            label72.Location = new System.Drawing.Point(302, 59);
            label72.Name = "label72";
            label72.Size = new System.Drawing.Size(109, 13);
            label72.TabIndex = 47;
            label72.Text = "fine clay translocation";
            // 
            // label68
            // 
            label68.AutoSize = true;
            label68.Location = new System.Drawing.Point(133, 103);
            label68.Name = "label68";
            label68.Size = new System.Drawing.Size(0, 13);
            label68.TabIndex = 59;
            // 
            // label71
            // 
            label71.AutoSize = true;
            label71.Location = new System.Drawing.Point(133, 77);
            label71.Name = "label71";
            label71.Size = new System.Drawing.Size(99, 13);
            label71.TabIndex = 58;
            label71.Text = "depth decay rate [-]";
            // 
            // label73
            // 
            label73.AutoSize = true;
            label73.Location = new System.Drawing.Point(133, 48);
            label73.Name = "label73";
            label73.Size = new System.Drawing.Size(167, 13);
            label73.TabIndex = 57;
            label73.Text = "potential bioturbation [kg / m2 / y]";
            // 
            // label74
            // 
            label74.AutoSize = true;
            label74.Location = new System.Drawing.Point(130, 117);
            label74.Name = "label74";
            label74.Size = new System.Drawing.Size(0, 13);
            label74.TabIndex = 64;
            // 
            // label75
            // 
            label75.AutoSize = true;
            label75.Location = new System.Drawing.Point(130, 91);
            label75.Name = "label75";
            label75.Size = new System.Drawing.Size(124, 13);
            label75.TabIndex = 63;
            label75.Text = "depth limitation rate [m-1]";
            // 
            // label76
            // 
            label76.AutoSize = true;
            label76.Location = new System.Drawing.Point(130, 62);
            label76.Name = "label76";
            label76.Size = new System.Drawing.Size(205, 13);
            label76.TabIndex = 62;
            label76.Text = "potential organic matter input [kg / m2 / y]";
            // 
            // label77
            // 
            label77.AutoSize = true;
            label77.Location = new System.Drawing.Point(130, 172);
            label77.Name = "label77";
            label77.Size = new System.Drawing.Size(0, 13);
            label77.TabIndex = 69;
            // 
            // label81
            // 
            label81.AutoSize = true;
            label81.Location = new System.Drawing.Point(130, 117);
            label81.Name = "label81";
            label81.Size = new System.Drawing.Size(113, 13);
            label81.TabIndex = 67;
            label81.Text = "humification fraction [-]";
            // 
            // label80
            // 
            label80.AutoSize = true;
            label80.Location = new System.Drawing.Point(381, 62);
            label80.Name = "label80";
            label80.Size = new System.Drawing.Size(133, 26);
            label80.TabIndex = 70;
            label80.Text = "decomposition parameters \r\nfor two OM pools:";
            // 
            // label82
            // 
            label82.AutoSize = true;
            label82.Location = new System.Drawing.Point(382, 94);
            label82.Name = "label82";
            label82.Size = new System.Drawing.Size(36, 13);
            label82.TabIndex = 71;
            label82.Text = "young";
            // 
            // label83
            // 
            label83.AutoSize = true;
            label83.Location = new System.Drawing.Point(496, 94);
            label83.Name = "label83";
            label83.Size = new System.Drawing.Size(21, 13);
            label83.TabIndex = 72;
            label83.Text = "old";
            // 
            // label84
            // 
            label84.AutoSize = true;
            label84.Location = new System.Drawing.Point(558, 140);
            label84.Name = "label84";
            label84.Size = new System.Drawing.Size(136, 13);
            label84.TabIndex = 74;
            label84.Text = "depth decay constant [m-1]";
            // 
            // label85
            // 
            label85.AutoSize = true;
            label85.Location = new System.Drawing.Point(558, 166);
            label85.Name = "label85";
            label85.Size = new System.Drawing.Size(116, 13);
            label85.TabIndex = 77;
            label85.Text = "TWI decay constant [-]";
            // 
            // label86
            // 
            label86.AutoSize = true;
            label86.Location = new System.Drawing.Point(558, 114);
            label86.Name = "label86";
            label86.Size = new System.Drawing.Size(115, 13);
            label86.TabIndex = 80;
            label86.Text = "decomposition rate [/y]";
            // 
            // label13
            // 
            label13.AutoSize = true;
            label13.Location = new System.Drawing.Point(410, 172);
            label13.Name = "label13";
            label13.Size = new System.Drawing.Size(112, 13);
            label13.TabIndex = 56;
            label13.Text = "Depth decay constant";
            label13.Click += new System.EventHandler(this.label13_Click);
            // 
            // mainMenu1
            // 
            this.mainMenu1.MenuItems.AddRange(new System.Windows.Forms.MenuItem[] {
            this.menuItemConfigFile,
            this.Mapselector,
            this.Menu_About_box});
            // 
            // menuItemConfigFile
            // 
            this.menuItemConfigFile.Index = 0;
            this.menuItemConfigFile.MenuItems.AddRange(new System.Windows.Forms.MenuItem[] {
            this.menuItemConfigFileOpen,
            this.menuItemConfigFileSaveAs,
            this.menuItemConfigFileSave});
            this.menuItemConfigFile.Text = "&RunFile";
            // 
            // menuItemConfigFileOpen
            // 
            this.menuItemConfigFileOpen.Index = 0;
            this.menuItemConfigFileOpen.Text = "&Open";
            this.menuItemConfigFileOpen.Click += new System.EventHandler(this.menuItemConfigFileOpen_Click);
            // 
            // menuItemConfigFileSaveAs
            // 
            this.menuItemConfigFileSaveAs.Index = 1;
            this.menuItemConfigFileSaveAs.Text = "Save &As";
            this.menuItemConfigFileSaveAs.Click += new System.EventHandler(this.menuItemConfigFileSave_Click);
            // 
            // menuItemConfigFileSave
            // 
            this.menuItemConfigFileSave.Index = 2;
            this.menuItemConfigFileSave.Text = "&Save";
            this.menuItemConfigFileSave.Click += new System.EventHandler(this.menuItemConfigFileSave_Click);
            // 
            // Mapselector
            // 
            this.Mapselector.Index = 1;
            this.Mapselector.MenuItems.AddRange(new System.Windows.Forms.MenuItem[] {
            this.Menu_map_total_sediment,
            this.Menu_map_waterflow,
            this.Menu_map_tillage,
            this.Menu_map_water_ero,
            this.Menu_map_creep,
            this.Menu_map_weathering,
            this.Menu_map_landsliding,
            this.Menu_map_critical_rainfall});
            this.Mapselector.Text = "&Map";
            // 
            // Menu_map_total_sediment
            // 
            this.Menu_map_total_sediment.Index = 0;
            this.Menu_map_total_sediment.Text = "Total redistribution";
            this.Menu_map_total_sediment.Click += new System.EventHandler(this.Menu_map_sediment_Click);
            // 
            // Menu_map_waterflow
            // 
            this.Menu_map_waterflow.Index = 1;
            this.Menu_map_waterflow.Text = "Annual water flow ";
            this.Menu_map_waterflow.Click += new System.EventHandler(this.Menu_map_waterflow_Click);
            // 
            // Menu_map_tillage
            // 
            this.Menu_map_tillage.Index = 2;
            this.Menu_map_tillage.Text = "Tillage";
            this.Menu_map_tillage.Click += new System.EventHandler(this.Menu_map_tillage_Click);
            // 
            // Menu_map_water_ero
            // 
            this.Menu_map_water_ero.Index = 3;
            this.Menu_map_water_ero.Text = "Overland erosion";
            this.Menu_map_water_ero.Click += new System.EventHandler(this.Menu_map_water_ero_Click);
            // 
            // Menu_map_creep
            // 
            this.Menu_map_creep.Index = 4;
            this.Menu_map_creep.Text = "Creep";
            this.Menu_map_creep.Click += new System.EventHandler(this.Menu_map_creep_Click);
            // 
            // Menu_map_weathering
            // 
            this.Menu_map_weathering.Index = 5;
            this.Menu_map_weathering.Text = "Weathering";
            this.Menu_map_weathering.Click += new System.EventHandler(this.Menu_map_landsliding_Click);
            // 
            // Menu_map_landsliding
            // 
            this.Menu_map_landsliding.Index = 6;
            this.Menu_map_landsliding.Text = "Landsliding";
            // 
            // Menu_map_critical_rainfall
            // 
            this.Menu_map_critical_rainfall.Index = 7;
            this.Menu_map_critical_rainfall.Text = "Cricital rainfall";
            this.Menu_map_critical_rainfall.Click += new System.EventHandler(this.Menu_map_critical_rainfall_Click);
            // 
            // Menu_About_box
            // 
            this.Menu_About_box.Index = 2;
            this.Menu_About_box.Text = "&About";
            this.Menu_About_box.Click += new System.EventHandler(this.Menu_aboutbox_Click);
            // 
            // statusBar1
            // 
            this.statusBar1.Location = new System.Drawing.Point(0, 475);
            this.statusBar1.Name = "statusBar1";
            this.statusBar1.Panels.AddRange(new System.Windows.Forms.StatusBarPanel[] {
            this.InfoStatusPanel,
            this.TimeStatusPanel,
            this.ProcessStatusPanel,
            this.out_sed_statuspanel,
            this.total_tillage_statuspanel});
            this.statusBar1.ShowPanels = true;
            this.statusBar1.Size = new System.Drawing.Size(1184, 22);
            this.statusBar1.SizingGrip = false;
            this.statusBar1.TabIndex = 144;
            this.statusBar1.Text = "statusBar1";
            this.statusBar1.PanelClick += new System.Windows.Forms.StatusBarPanelClickEventHandler(this.statusBar1_PanelClick);
            // 
            // InfoStatusPanel
            // 
            this.InfoStatusPanel.Name = "InfoStatusPanel";
            this.InfoStatusPanel.Text = "info";
            this.InfoStatusPanel.Width = 200;
            // 
            // TimeStatusPanel
            // 
            this.TimeStatusPanel.Name = "TimeStatusPanel";
            this.TimeStatusPanel.Text = "time";
            this.TimeStatusPanel.Width = 80;
            // 
            // ProcessStatusPanel
            // 
            this.ProcessStatusPanel.Name = "ProcessStatusPanel";
            this.ProcessStatusPanel.Text = "processes";
            this.ProcessStatusPanel.Width = 120;
            // 
            // out_sed_statuspanel
            // 
            this.out_sed_statuspanel.Name = "out_sed_statuspanel";
            this.out_sed_statuspanel.Text = "sed export";
            this.out_sed_statuspanel.Width = 140;
            // 
            // total_tillage_statuspanel
            // 
            this.total_tillage_statuspanel.Name = "total_tillage_statuspanel";
            this.total_tillage_statuspanel.Text = "tillage volume";
            this.total_tillage_statuspanel.Width = 140;
            // 
            // start_button
            // 
            this.start_button.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Left)));
            this.start_button.Location = new System.Drawing.Point(9, 419);
            this.start_button.Name = "start_button";
            this.start_button.Size = new System.Drawing.Size(88, 27);
            this.start_button.TabIndex = 146;
            this.start_button.Text = "Start";
            this.start_button.Click += new System.EventHandler(this.start_run);
            // 
            // End_button
            // 
            this.End_button.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Left)));
            this.End_button.Location = new System.Drawing.Point(103, 419);
            this.End_button.Name = "End_button";
            this.End_button.Size = new System.Drawing.Size(100, 27);
            this.End_button.TabIndex = 147;
            this.End_button.Text = "Quit";
            this.End_button.Click += new System.EventHandler(this.End_button_Click);
            // 
            // View_tabs_checkbox
            // 
            this.View_tabs_checkbox.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Left)));
            this.View_tabs_checkbox.AutoSize = true;
            this.View_tabs_checkbox.Checked = true;
            this.View_tabs_checkbox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.View_tabs_checkbox.Location = new System.Drawing.Point(354, 427);
            this.View_tabs_checkbox.Name = "View_tabs_checkbox";
            this.View_tabs_checkbox.Size = new System.Drawing.Size(77, 17);
            this.View_tabs_checkbox.TabIndex = 151;
            this.View_tabs_checkbox.Text = "view tabs?";
            this.View_tabs_checkbox.UseVisualStyleBackColor = true;
            this.View_tabs_checkbox.CheckedChanged += new System.EventHandler(this.View_tabs_checkbox_CheckedChanged);
            // 
            // label2
            // 
            this.label2.Location = new System.Drawing.Point(19, 121);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(104, 24);
            this.label2.TabIndex = 97;
            this.label2.Text = "Rainfall data file";
            this.label2.TextAlign = System.Drawing.ContentAlignment.MiddleRight;
            this.toolTip1.SetToolTip(this.label2, "Hourly rainfall data - in an ascii format");
            // 
            // checkBoxGenerateAVIFile
            // 
            this.checkBoxGenerateAVIFile.Location = new System.Drawing.Point(18, 26);
            this.checkBoxGenerateAVIFile.Name = "checkBoxGenerateAVIFile";
            this.checkBoxGenerateAVIFile.Size = new System.Drawing.Size(128, 24);
            this.checkBoxGenerateAVIFile.TabIndex = 200;
            this.checkBoxGenerateAVIFile.Text = "Generate Avi File";
            this.toolTip1.SetToolTip(this.checkBoxGenerateAVIFile, "Check to generate a movie file of the screen display");
            // 
            // label33
            // 
            this.label33.Location = new System.Drawing.Point(18, 53);
            this.label33.Name = "label33";
            this.label33.Size = new System.Drawing.Size(131, 25);
            this.label33.TabIndex = 203;
            this.label33.Text = "Save map every * years";
            this.label33.TextAlign = System.Drawing.ContentAlignment.MiddleRight;
            this.toolTip1.SetToolTip(this.label33, "How often the avi file AND the other data files are saved");
            // 
            // googleAnimationCheckbox
            // 
            this.googleAnimationCheckbox.Location = new System.Drawing.Point(18, 86);
            this.googleAnimationCheckbox.Name = "googleAnimationCheckbox";
            this.googleAnimationCheckbox.Size = new System.Drawing.Size(195, 24);
            this.googleAnimationCheckbox.TabIndex = 210;
            this.googleAnimationCheckbox.Text = "Google Earth Animation ";
            this.toolTip1.SetToolTip(this.googleAnimationCheckbox, "Check to generate a movie file of the screen display");
            // 
            // label79
            // 
            this.label79.Location = new System.Drawing.Point(15, 111);
            this.label79.Name = "label79";
            this.label79.Size = new System.Drawing.Size(134, 25);
            this.label79.TabIndex = 213;
            this.label79.Text = "Save map every * years";
            this.label79.TextAlign = System.Drawing.ContentAlignment.MiddleRight;
            this.toolTip1.SetToolTip(this.label79, "How often the avi file AND the other data files are saved");
            // 
            // textBoxAVIFile
            // 
            this.textBoxAVIFile.Location = new System.Drawing.Point(165, 29);
            this.textBoxAVIFile.Name = "textBoxAVIFile";
            this.textBoxAVIFile.Size = new System.Drawing.Size(112, 20);
            this.textBoxAVIFile.TabIndex = 199;
            this.textBoxAVIFile.Text = "out.avi";
            this.toolTip1.SetToolTip(this.textBoxAVIFile, "File name for avi file");
            // 
            // googleAnimationTextBox
            // 
            this.googleAnimationTextBox.Location = new System.Drawing.Point(165, 88);
            this.googleAnimationTextBox.Name = "googleAnimationTextBox";
            this.googleAnimationTextBox.Size = new System.Drawing.Size(112, 20);
            this.googleAnimationTextBox.TabIndex = 211;
            this.googleAnimationTextBox.Text = "animation.kmz";
            this.toolTip1.SetToolTip(this.googleAnimationTextBox, "File name for avi file");
            // 
            // UTMzonebox
            // 
            this.UTMzonebox.Location = new System.Drawing.Point(99, 16);
            this.UTMzonebox.Name = "UTMzonebox";
            this.UTMzonebox.Size = new System.Drawing.Size(39, 20);
            this.UTMzonebox.TabIndex = 194;
            this.UTMzonebox.Tag = "UTM zone";
            this.toolTip1.SetToolTip(this.UTMzonebox, "Enter the UTM zone");
            // 
            // label25
            // 
            this.label25.Location = new System.Drawing.Point(19, 164);
            this.label25.Name = "label25";
            this.label25.Size = new System.Drawing.Size(143, 24);
            this.label25.TabIndex = 97;
            this.label25.Text = "mean annual rainfall [m]";
            this.label25.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            this.toolTip1.SetToolTip(this.label25, "Hourly rainfall data - in an ascii format");
            // 
            // label14
            // 
            this.label14.Location = new System.Drawing.Point(19, 190);
            this.label14.Name = "label14";
            this.label14.Size = new System.Drawing.Size(143, 24);
            this.label14.TabIndex = 115;
            this.label14.Text = "mean annual infiltration [m]";
            this.label14.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            this.toolTip1.SetToolTip(this.label14, "Hourly rainfall data - in an ascii format");
            // 
            // label15
            // 
            this.label15.Location = new System.Drawing.Point(19, 216);
            this.label15.Name = "label15";
            this.label15.Size = new System.Drawing.Size(151, 24);
            this.label15.TabIndex = 118;
            this.label15.Text = "mean annual evaporation [m]";
            this.label15.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            this.toolTip1.SetToolTip(this.label15, "Hourly rainfall data - in an ascii format");
            // 
            // label17
            // 
            this.label17.Location = new System.Drawing.Point(19, 140);
            this.label17.Name = "label17";
            this.label17.Size = new System.Drawing.Size(143, 24);
            this.label17.TabIndex = 121;
            this.label17.Text = "tillage fields [1/-9999]";
            this.label17.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            this.toolTip1.SetToolTip(this.label17, "Hourly rainfall data - in an ascii format");
            // 
            // fill_sinks_before_checkbox
            // 
            this.fill_sinks_before_checkbox.AutoSize = true;
            this.fill_sinks_before_checkbox.Location = new System.Drawing.Point(11, 22);
            this.fill_sinks_before_checkbox.Name = "fill_sinks_before_checkbox";
            this.fill_sinks_before_checkbox.Size = new System.Drawing.Size(131, 17);
            this.fill_sinks_before_checkbox.TabIndex = 132;
            this.fill_sinks_before_checkbox.Text = "remove sinks and flats";
            this.toolTip1.SetToolTip(this.fill_sinks_before_checkbox, resources.GetString("fill_sinks_before_checkbox.ToolTip"));
            this.fill_sinks_before_checkbox.UseVisualStyleBackColor = true;
            // 
            // label8
            // 
            this.label8.Location = new System.Drawing.Point(110, 31);
            this.label8.Name = "label8";
            this.label8.Size = new System.Drawing.Size(54, 25);
            this.label8.TabIndex = 222;
            this.label8.Text = "timesteps";
            this.label8.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            this.toolTip1.SetToolTip(this.label8, "How often the avi file AND the other data files are saved");
            // 
            // fill_sinks_during_checkbox
            // 
            this.fill_sinks_during_checkbox.AutoSize = true;
            this.fill_sinks_during_checkbox.Location = new System.Drawing.Point(11, 25);
            this.fill_sinks_during_checkbox.Name = "fill_sinks_during_checkbox";
            this.fill_sinks_during_checkbox.Size = new System.Drawing.Size(131, 17);
            this.fill_sinks_during_checkbox.TabIndex = 132;
            this.fill_sinks_during_checkbox.Text = "remove sinks and flats";
            this.toolTip1.SetToolTip(this.fill_sinks_during_checkbox, resources.GetString("fill_sinks_during_checkbox.ToolTip"));
            this.fill_sinks_during_checkbox.UseVisualStyleBackColor = true;
            // 
            // groupBox13
            // 
            this.groupBox13.Controls.Add(this.fill_sinks_during_checkbox);
            this.groupBox13.Location = new System.Drawing.Point(589, 91);
            this.groupBox13.Name = "groupBox13";
            this.groupBox13.Size = new System.Drawing.Size(158, 55);
            this.groupBox13.TabIndex = 135;
            this.groupBox13.TabStop = false;
            this.groupBox13.Text = "while running: ";
            this.toolTip1.SetToolTip(this.groupBox13, "LORICA");
            // 
            // label7
            // 
            this.label7.Location = new System.Drawing.Point(407, 22);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(134, 24);
            this.label7.TabIndex = 110;
            this.label7.Text = " constant value";
            this.label7.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            this.toolTip1.SetToolTip(this.label7, "\"constant value\" is only available for some inputs and only \r\nwhen neither f(t) n" +
        "or f(x,y) are checked");
            // 
            // label5
            // 
            this.label5.Location = new System.Drawing.Point(258, 23);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(51, 24);
            this.label5.TabIndex = 108;
            this.label5.Text = "filename";
            this.label5.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            this.toolTip1.SetToolTip(this.label5, resources.GetString("label5.ToolTip"));
            // 
            // label29
            // 
            this.label29.Location = new System.Drawing.Point(217, 22);
            this.label29.Name = "label29";
            this.label29.Size = new System.Drawing.Size(24, 24);
            this.label29.TabIndex = 139;
            this.label29.Text = "f(t)";
            this.label29.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            this.toolTip1.SetToolTip(this.label29, "Check this for temporally variable inputs");
            // 
            // groupBox3
            // 
            this.groupBox3.Controls.Add(this.landuse_determinator_button);
            this.groupBox3.Location = new System.Drawing.Point(589, 152);
            this.groupBox3.Name = "groupBox3";
            this.groupBox3.Size = new System.Drawing.Size(158, 55);
            this.groupBox3.TabIndex = 136;
            this.groupBox3.TabStop = false;
            this.groupBox3.Text = "for landuse: ";
            this.toolTip1.SetToolTip(this.groupBox3, "LORICA");
            // 
            // landuse_determinator_button
            // 
            this.landuse_determinator_button.Location = new System.Drawing.Point(11, 19);
            this.landuse_determinator_button.Name = "landuse_determinator_button";
            this.landuse_determinator_button.Size = new System.Drawing.Size(122, 23);
            this.landuse_determinator_button.TabIndex = 0;
            this.landuse_determinator_button.Text = "determine effects";
            this.landuse_determinator_button.UseVisualStyleBackColor = true;
            this.landuse_determinator_button.Click += new System.EventHandler(this.landuse_determinator_button_Click);
            // 
            // groupBox9
            // 
            this.groupBox9.Controls.Add(this.parameter_k1_textbox);
            this.groupBox9.Controls.Add(this.label24);
            this.groupBox9.Controls.Add(this.label26);
            this.groupBox9.Controls.Add(this.label27);
            this.groupBox9.Controls.Add(this.label28);
            this.groupBox9.Controls.Add(this.parameter_k2_textbox);
            this.groupBox9.Controls.Add(this.parameter_Pa_textbox);
            this.groupBox9.Controls.Add(this.parameter_P0_textbox);
            this.groupBox9.Controls.Add(this.label21);
            this.groupBox9.Controls.Add(this.Biological_weathering_checkbox);
            this.groupBox9.Location = new System.Drawing.Point(12, 13);
            this.groupBox9.Name = "groupBox9";
            this.groupBox9.Size = new System.Drawing.Size(222, 180);
            this.groupBox9.TabIndex = 5;
            this.groupBox9.TabStop = false;
            this.groupBox9.Text = "Biological weathering (humped model)";
            this.toolTip1.SetToolTip(this.groupBox9, resources.GetString("groupBox9.ToolTip"));
            // 
            // parameter_k1_textbox
            // 
            this.parameter_k1_textbox.Location = new System.Drawing.Point(14, 80);
            this.parameter_k1_textbox.Name = "parameter_k1_textbox";
            this.parameter_k1_textbox.Size = new System.Drawing.Size(53, 20);
            this.parameter_k1_textbox.TabIndex = 20;
            this.parameter_k1_textbox.Text = "0.1";
            // 
            // label24
            // 
            this.label24.AutoSize = true;
            this.label24.Location = new System.Drawing.Point(89, 135);
            this.label24.Name = "label24";
            this.label24.Size = new System.Drawing.Size(52, 13);
            this.label24.TabIndex = 19;
            this.label24.Text = "Pa (m t-1)";
            // 
            // label26
            // 
            this.label26.AutoSize = true;
            this.label26.Location = new System.Drawing.Point(89, 109);
            this.label26.Name = "label26";
            this.label26.Size = new System.Drawing.Size(40, 13);
            this.label26.TabIndex = 18;
            this.label26.Text = "k2 (t-1)";
            // 
            // label27
            // 
            this.label27.AutoSize = true;
            this.label27.Location = new System.Drawing.Point(89, 83);
            this.label27.Name = "label27";
            this.label27.Size = new System.Drawing.Size(40, 13);
            this.label27.TabIndex = 17;
            this.label27.Text = "k1 (t-1)";
            // 
            // label28
            // 
            this.label28.AutoSize = true;
            this.label28.Location = new System.Drawing.Point(89, 57);
            this.label28.Name = "label28";
            this.label28.Size = new System.Drawing.Size(52, 13);
            this.label28.TabIndex = 16;
            this.label28.Text = "P0 (m t-1)";
            // 
            // parameter_k2_textbox
            // 
            this.parameter_k2_textbox.Location = new System.Drawing.Point(14, 106);
            this.parameter_k2_textbox.Name = "parameter_k2_textbox";
            this.parameter_k2_textbox.Size = new System.Drawing.Size(53, 20);
            this.parameter_k2_textbox.TabIndex = 14;
            this.parameter_k2_textbox.Text = "6";
            // 
            // parameter_Pa_textbox
            // 
            this.parameter_Pa_textbox.Location = new System.Drawing.Point(14, 132);
            this.parameter_Pa_textbox.Name = "parameter_Pa_textbox";
            this.parameter_Pa_textbox.Size = new System.Drawing.Size(53, 20);
            this.parameter_Pa_textbox.TabIndex = 13;
            this.parameter_Pa_textbox.Text = "0.00002";
            // 
            // parameter_P0_textbox
            // 
            this.parameter_P0_textbox.Location = new System.Drawing.Point(14, 54);
            this.parameter_P0_textbox.Name = "parameter_P0_textbox";
            this.parameter_P0_textbox.Size = new System.Drawing.Size(53, 20);
            this.parameter_P0_textbox.TabIndex = 12;
            this.parameter_P0_textbox.Text = "0.000033";
            // 
            // label21
            // 
            this.label21.AutoSize = true;
            this.label21.Location = new System.Drawing.Point(11, 11);
            this.label21.Name = "label21";
            this.label21.Size = new System.Drawing.Size(0, 13);
            this.label21.TabIndex = 4;
            // 
            // Biological_weathering_checkbox
            // 
            this.Biological_weathering_checkbox.AutoSize = true;
            this.Biological_weathering_checkbox.Location = new System.Drawing.Point(14, 19);
            this.Biological_weathering_checkbox.Name = "Biological_weathering_checkbox";
            this.Biological_weathering_checkbox.Size = new System.Drawing.Size(124, 17);
            this.Biological_weathering_checkbox.TabIndex = 3;
            this.Biological_weathering_checkbox.Text = "Activate this process";
            this.Biological_weathering_checkbox.UseVisualStyleBackColor = true;
            // 
            // label98
            // 
            this.label98.Location = new System.Drawing.Point(19, 245);
            this.label98.Name = "label98";
            this.label98.Size = new System.Drawing.Size(151, 24);
            this.label98.TabIndex = 152;
            this.label98.Text = "mean annual temperature [C]";
            this.label98.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            this.toolTip1.SetToolTip(this.label98, "Hourly rainfall data - in an ascii format");
            this.label98.Click += new System.EventHandler(this.label98_Click);
            // 
            // label99
            // 
            this.label99.Location = new System.Drawing.Point(52, 269);
            this.label99.Name = "label99";
            this.label99.Size = new System.Drawing.Size(151, 24);
            this.label99.TabIndex = 153;
            this.label99.Text = "Only with daily water balance";
            this.label99.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            this.toolTip1.SetToolTip(this.label99, "Hourly rainfall data - in an ascii format");
            // 
            // trackBar1
            // 
            this.trackBar1.AutoSize = false;
            this.trackBar1.Enabled = false;
            this.trackBar1.Location = new System.Drawing.Point(215, 19);
            this.trackBar1.Name = "trackBar1";
            this.trackBar1.Size = new System.Drawing.Size(104, 20);
            this.trackBar1.TabIndex = 149;
            this.trackBar1.Scroll += new System.EventHandler(this.trackBar1_Scroll);
            // 
            // label61
            // 
            this.label61.AutoSize = true;
            this.label61.Location = new System.Drawing.Point(5, 24);
            this.label61.Name = "label61";
            this.label61.Size = new System.Drawing.Size(44, 13);
            this.label61.TabIndex = 150;
            this.label61.Text = "Graphic";
            this.label61.Visible = false;
            // 
            // map_controls
            // 
            this.map_controls.Controls.Add(this.label63);
            this.map_controls.Controls.Add(this.graphicToGoogleEarthButton);
            this.map_controls.Controls.Add(this.trackBar2);
            this.map_controls.Controls.Add(this.label62);
            this.map_controls.Controls.Add(this.comboBox1);
            this.map_controls.Controls.Add(this.trackBar1);
            this.map_controls.Controls.Add(this.label61);
            this.map_controls.Location = new System.Drawing.Point(481, 372);
            this.map_controls.Name = "map_controls";
            this.map_controls.Size = new System.Drawing.Size(334, 91);
            this.map_controls.TabIndex = 151;
            this.map_controls.TabStop = false;
            this.map_controls.Text = "Map controls";
            this.map_controls.Visible = false;
            // 
            // label63
            // 
            this.label63.AutoSize = true;
            this.label63.Location = new System.Drawing.Point(175, 55);
            this.label63.Name = "label63";
            this.label63.Size = new System.Drawing.Size(34, 13);
            this.label63.TabIndex = 154;
            this.label63.Text = "Zoom";
            // 
            // graphicToGoogleEarthButton
            // 
            this.graphicToGoogleEarthButton.Location = new System.Drawing.Point(26, 49);
            this.graphicToGoogleEarthButton.Name = "graphicToGoogleEarthButton";
            this.graphicToGoogleEarthButton.Size = new System.Drawing.Size(131, 25);
            this.graphicToGoogleEarthButton.TabIndex = 150;
            this.graphicToGoogleEarthButton.Text = "image to google earth";
            this.graphicToGoogleEarthButton.UseVisualStyleBackColor = true;
            this.graphicToGoogleEarthButton.Visible = false;
            // 
            // trackBar2
            // 
            this.trackBar2.AutoSize = false;
            this.trackBar2.Enabled = false;
            this.trackBar2.Location = new System.Drawing.Point(215, 54);
            this.trackBar2.Minimum = 1;
            this.trackBar2.Name = "trackBar2";
            this.trackBar2.Size = new System.Drawing.Size(104, 20);
            this.trackBar2.TabIndex = 153;
            this.trackBar2.Value = 5;
            this.trackBar2.Scroll += new System.EventHandler(this.trackBar2_Scroll);
            // 
            // label62
            // 
            this.label62.AutoSize = true;
            this.label62.Location = new System.Drawing.Point(169, 24);
            this.label62.Name = "label62";
            this.label62.Size = new System.Drawing.Size(46, 13);
            this.label62.TabIndex = 152;
            this.label62.Text = "Contrast";
            // 
            // comboBox1
            // 
            this.comboBox1.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            this.comboBox1.FormattingEnabled = true;
            this.comboBox1.Location = new System.Drawing.Point(53, 20);
            this.comboBox1.Name = "comboBox1";
            this.comboBox1.Size = new System.Drawing.Size(104, 21);
            this.comboBox1.TabIndex = 151;
            this.comboBox1.Visible = false;
            // 
            // openFileDialog1
            // 
            this.openFileDialog1.FileName = "openFileDialog1";
            // 
            // label1
            // 
            this.label1.Location = new System.Drawing.Point(19, 49);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(104, 24);
            this.label1.TabIndex = 56;
            this.label1.Text = "DEM (.asc format)";
            this.label1.TextAlign = System.Drawing.ContentAlignment.MiddleRight;
            // 
            // button6
            // 
            this.button6.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Left)));
            this.button6.Location = new System.Drawing.Point(132, 235);
            this.button6.Name = "button6";
            this.button6.Size = new System.Drawing.Size(100, 24);
            this.button6.TabIndex = 7;
            this.button6.Text = "test data";
            // 
            // textBox1
            // 
            this.textBox1.Location = new System.Drawing.Point(131, 121);
            this.textBox1.Name = "textBox1";
            this.textBox1.Size = new System.Drawing.Size(120, 20);
            this.textBox1.TabIndex = 103;
            this.textBox1.Text = "null";
            // 
            // textBox2
            // 
            this.textBox2.Location = new System.Drawing.Point(131, 49);
            this.textBox2.Name = "textBox2";
            this.textBox2.Size = new System.Drawing.Size(120, 20);
            this.textBox2.TabIndex = 100;
            this.textBox2.Text = "whole9.dat";
            // 
            // Output
            // 
            this.Output.Controls.Add(this.groupBox6);
            this.Output.Controls.Add(this.groupBox5);
            this.Output.Location = new System.Drawing.Point(4, 22);
            this.Output.Name = "Output";
            this.Output.Size = new System.Drawing.Size(803, 293);
            this.Output.TabIndex = 7;
            this.Output.Text = "Output";
            this.Output.UseVisualStyleBackColor = true;
            // 
            // groupBox6
            // 
            this.groupBox6.Controls.Add(this.groupBox12);
            this.groupBox6.Controls.Add(this.groupBox11);
            this.groupBox6.Controls.Add(this.groupBox1);
            this.groupBox6.Controls.Add(this.checkedListBox1);
            this.groupBox6.Location = new System.Drawing.Point(8, 16);
            this.groupBox6.Name = "groupBox6";
            this.groupBox6.Size = new System.Drawing.Size(294, 258);
            this.groupBox6.TabIndex = 222;
            this.groupBox6.TabStop = false;
            this.groupBox6.Text = "Normal outputs (ascii grids)";
            // 
            // groupBox12
            // 
            this.groupBox12.Controls.Add(this.annual_output_checkbox);
            this.groupBox12.Controls.Add(this.cumulative_output_checkbox);
            this.groupBox12.Location = new System.Drawing.Point(181, 19);
            this.groupBox12.Name = "groupBox12";
            this.groupBox12.Size = new System.Drawing.Size(102, 63);
            this.groupBox12.TabIndex = 227;
            this.groupBox12.TabStop = false;
            // 
            // annual_output_checkbox
            // 
            this.annual_output_checkbox.AutoSize = true;
            this.annual_output_checkbox.Location = new System.Drawing.Point(5, 35);
            this.annual_output_checkbox.Name = "annual_output_checkbox";
            this.annual_output_checkbox.Size = new System.Drawing.Size(57, 17);
            this.annual_output_checkbox.TabIndex = 1;
            this.annual_output_checkbox.Text = "annual";
            this.annual_output_checkbox.UseVisualStyleBackColor = true;
            // 
            // cumulative_output_checkbox
            // 
            this.cumulative_output_checkbox.AutoSize = true;
            this.cumulative_output_checkbox.Checked = true;
            this.cumulative_output_checkbox.Location = new System.Drawing.Point(5, 12);
            this.cumulative_output_checkbox.Name = "cumulative_output_checkbox";
            this.cumulative_output_checkbox.Size = new System.Drawing.Size(76, 17);
            this.cumulative_output_checkbox.TabIndex = 0;
            this.cumulative_output_checkbox.TabStop = true;
            this.cumulative_output_checkbox.Text = "cumulative";
            this.cumulative_output_checkbox.UseVisualStyleBackColor = true;
            // 
            // groupBox11
            // 
            this.groupBox11.Controls.Add(this.label8);
            this.groupBox11.Controls.Add(this.Regular_output_checkbox);
            this.groupBox11.Controls.Add(this.Final_output_checkbox);
            this.groupBox11.Controls.Add(this.Box_years_output);
            this.groupBox11.Location = new System.Drawing.Point(6, 19);
            this.groupBox11.Name = "groupBox11";
            this.groupBox11.Size = new System.Drawing.Size(166, 63);
            this.groupBox11.TabIndex = 226;
            this.groupBox11.TabStop = false;
            // 
            // Regular_output_checkbox
            // 
            this.Regular_output_checkbox.AutoSize = true;
            this.Regular_output_checkbox.Location = new System.Drawing.Point(6, 36);
            this.Regular_output_checkbox.Name = "Regular_output_checkbox";
            this.Regular_output_checkbox.Size = new System.Drawing.Size(55, 17);
            this.Regular_output_checkbox.TabIndex = 221;
            this.Regular_output_checkbox.Text = "every ";
            this.Regular_output_checkbox.UseVisualStyleBackColor = true;
            // 
            // Final_output_checkbox
            // 
            this.Final_output_checkbox.AutoSize = true;
            this.Final_output_checkbox.Checked = true;
            this.Final_output_checkbox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.Final_output_checkbox.Location = new System.Drawing.Point(6, 13);
            this.Final_output_checkbox.Name = "Final_output_checkbox";
            this.Final_output_checkbox.Size = new System.Drawing.Size(81, 17);
            this.Final_output_checkbox.TabIndex = 220;
            this.Final_output_checkbox.Text = "when ready";
            this.Final_output_checkbox.UseVisualStyleBackColor = true;
            // 
            // Box_years_output
            // 
            this.Box_years_output.AcceptsTab = true;
            this.Box_years_output.Location = new System.Drawing.Point(67, 34);
            this.Box_years_output.Name = "Box_years_output";
            this.Box_years_output.Size = new System.Drawing.Size(44, 20);
            this.Box_years_output.TabIndex = 1;
            this.Box_years_output.Text = "3";
            // 
            // groupBox1
            // 
            this.groupBox1.Controls.Add(this.diagnostic_output_checkbox);
            this.groupBox1.Controls.Add(this.label37);
            this.groupBox1.Controls.Add(this.outputcode_textbox);
            this.groupBox1.Controls.Add(this.water_output_checkbox);
            this.groupBox1.Controls.Add(this.depressions_output_checkbox);
            this.groupBox1.Controls.Add(this.all_process_output_checkbox);
            this.groupBox1.Controls.Add(this.Soildepth_output_checkbox);
            this.groupBox1.Controls.Add(this.Alt_change_output_checkbox);
            this.groupBox1.Controls.Add(this.Altitude_output_checkbox);
            this.groupBox1.Location = new System.Drawing.Point(6, 85);
            this.groupBox1.Name = "groupBox1";
            this.groupBox1.Size = new System.Drawing.Size(277, 167);
            this.groupBox1.TabIndex = 225;
            this.groupBox1.TabStop = false;
            this.groupBox1.Text = "Outputs:";
            // 
            // diagnostic_output_checkbox
            // 
            this.diagnostic_output_checkbox.AutoSize = true;
            this.diagnostic_output_checkbox.Font = new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Italic, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.diagnostic_output_checkbox.Location = new System.Drawing.Point(125, 141);
            this.diagnostic_output_checkbox.Name = "diagnostic_output_checkbox";
            this.diagnostic_output_checkbox.Size = new System.Drawing.Size(81, 17);
            this.diagnostic_output_checkbox.TabIndex = 230;
            this.diagnostic_output_checkbox.Text = "Diagnostics";
            this.diagnostic_output_checkbox.UseVisualStyleBackColor = true;
            // 
            // label37
            // 
            this.label37.AutoSize = true;
            this.label37.Location = new System.Drawing.Point(153, 31);
            this.label37.Name = "label37";
            this.label37.Size = new System.Drawing.Size(69, 13);
            this.label37.TabIndex = 229;
            this.label37.Text = "Output code:";
            // 
            // outputcode_textbox
            // 
            this.outputcode_textbox.Location = new System.Drawing.Point(156, 47);
            this.outputcode_textbox.Name = "outputcode_textbox";
            this.outputcode_textbox.Size = new System.Drawing.Size(100, 20);
            this.outputcode_textbox.TabIndex = 228;
            // 
            // water_output_checkbox
            // 
            this.water_output_checkbox.AutoSize = true;
            this.water_output_checkbox.Checked = true;
            this.water_output_checkbox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.water_output_checkbox.Font = new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.water_output_checkbox.Location = new System.Drawing.Point(24, 118);
            this.water_output_checkbox.Name = "water_output_checkbox";
            this.water_output_checkbox.Size = new System.Drawing.Size(74, 17);
            this.water_output_checkbox.TabIndex = 227;
            this.water_output_checkbox.Text = "Waterflow";
            this.water_output_checkbox.UseVisualStyleBackColor = true;
            // 
            // depressions_output_checkbox
            // 
            this.depressions_output_checkbox.AutoSize = true;
            this.depressions_output_checkbox.Font = new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Italic, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.depressions_output_checkbox.Location = new System.Drawing.Point(24, 141);
            this.depressions_output_checkbox.Name = "depressions_output_checkbox";
            this.depressions_output_checkbox.Size = new System.Drawing.Size(84, 17);
            this.depressions_output_checkbox.TabIndex = 226;
            this.depressions_output_checkbox.Text = "Depressions";
            this.depressions_output_checkbox.UseVisualStyleBackColor = true;
            // 
            // all_process_output_checkbox
            // 
            this.all_process_output_checkbox.AutoSize = true;
            this.all_process_output_checkbox.Checked = true;
            this.all_process_output_checkbox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.all_process_output_checkbox.Location = new System.Drawing.Point(24, 95);
            this.all_process_output_checkbox.Name = "all_process_output_checkbox";
            this.all_process_output_checkbox.Size = new System.Drawing.Size(134, 17);
            this.all_process_output_checkbox.TabIndex = 225;
            this.all_process_output_checkbox.Text = "Indiv. process volumes";
            this.all_process_output_checkbox.UseVisualStyleBackColor = true;
            // 
            // Soildepth_output_checkbox
            // 
            this.Soildepth_output_checkbox.AutoSize = true;
            this.Soildepth_output_checkbox.Checked = true;
            this.Soildepth_output_checkbox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.Soildepth_output_checkbox.Location = new System.Drawing.Point(24, 72);
            this.Soildepth_output_checkbox.Name = "Soildepth_output_checkbox";
            this.Soildepth_output_checkbox.Size = new System.Drawing.Size(70, 17);
            this.Soildepth_output_checkbox.TabIndex = 224;
            this.Soildepth_output_checkbox.Text = "Soildepth";
            this.Soildepth_output_checkbox.UseVisualStyleBackColor = true;
            // 
            // Alt_change_output_checkbox
            // 
            this.Alt_change_output_checkbox.AutoSize = true;
            this.Alt_change_output_checkbox.Checked = true;
            this.Alt_change_output_checkbox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.Alt_change_output_checkbox.Location = new System.Drawing.Point(24, 49);
            this.Alt_change_output_checkbox.Name = "Alt_change_output_checkbox";
            this.Alt_change_output_checkbox.Size = new System.Drawing.Size(100, 17);
            this.Alt_change_output_checkbox.TabIndex = 223;
            this.Alt_change_output_checkbox.Text = "Altitude change";
            this.Alt_change_output_checkbox.UseVisualStyleBackColor = true;
            // 
            // Altitude_output_checkbox
            // 
            this.Altitude_output_checkbox.AutoSize = true;
            this.Altitude_output_checkbox.Checked = true;
            this.Altitude_output_checkbox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.Altitude_output_checkbox.Location = new System.Drawing.Point(24, 26);
            this.Altitude_output_checkbox.Name = "Altitude_output_checkbox";
            this.Altitude_output_checkbox.Size = new System.Drawing.Size(61, 17);
            this.Altitude_output_checkbox.TabIndex = 222;
            this.Altitude_output_checkbox.Text = "Altitude";
            this.Altitude_output_checkbox.UseVisualStyleBackColor = true;
            // 
            // checkedListBox1
            // 
            this.checkedListBox1.FormattingEnabled = true;
            this.checkedListBox1.Items.AddRange(new object[] {
            "altitude",
            "altitude change",
            "soildepth",
            "soildepth change",
            "redistribution all processes",
            "redistribution per process",
            "weathering all processes",
            "weathering per process"});
            this.checkedListBox1.Location = new System.Drawing.Point(119, 96);
            this.checkedListBox1.Name = "checkedListBox1";
            this.checkedListBox1.Size = new System.Drawing.Size(152, 49);
            this.checkedListBox1.TabIndex = 0;
            this.checkedListBox1.Visible = false;
            // 
            // groupBox5
            // 
            this.groupBox5.Controls.Add(this.button3);
            this.groupBox5.Controls.Add(this.button1);
            this.groupBox5.Controls.Add(this.button2);
            this.groupBox5.Controls.Add(this.timeseries_form_button);
            this.groupBox5.Controls.Add(this.UTMgroupBox);
            this.groupBox5.Controls.Add(this.UTMgridcheckbox);
            this.groupBox5.Controls.Add(this.textBox4);
            this.groupBox5.Controls.Add(this.googleBeginDate);
            this.groupBox5.Controls.Add(this.googAnimationSaveInterval);
            this.groupBox5.Controls.Add(this.googleAnimationTextBox);
            this.groupBox5.Controls.Add(this.saveintervalbox);
            this.groupBox5.Controls.Add(this.textBoxAVIFile);
            this.groupBox5.Controls.Add(this.label78);
            this.groupBox5.Controls.Add(this.label79);
            this.groupBox5.Controls.Add(this.googleAnimationCheckbox);
            this.groupBox5.Controls.Add(this.label33);
            this.groupBox5.Controls.Add(this.checkBoxGenerateAVIFile);
            this.groupBox5.Location = new System.Drawing.Point(308, 15);
            this.groupBox5.Name = "groupBox5";
            this.groupBox5.Size = new System.Drawing.Size(464, 259);
            this.groupBox5.TabIndex = 221;
            this.groupBox5.TabStop = false;
            this.groupBox5.Text = "Other outputs";
            // 
            // button3
            // 
            this.button3.Location = new System.Drawing.Point(175, 218);
            this.button3.Name = "button3";
            this.button3.Size = new System.Drawing.Size(142, 24);
            this.button3.TabIndex = 222;
            this.button3.Text = "Google Earth output..";
            this.button3.UseVisualStyleBackColor = true;
            this.button3.Visible = false;
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(175, 189);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(142, 24);
            this.button1.TabIndex = 221;
            this.button1.Text = "video output..";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Visible = false;
            // 
            // button2
            // 
            this.button2.Location = new System.Drawing.Point(18, 218);
            this.button2.Name = "button2";
            this.button2.Size = new System.Drawing.Size(142, 24);
            this.button2.TabIndex = 220;
            this.button2.Text = "profile outputs..";
            this.button2.UseVisualStyleBackColor = true;
            this.button2.Click += new System.EventHandler(this.profiles_button_Click);
            // 
            // timeseries_form_button
            // 
            this.timeseries_form_button.Location = new System.Drawing.Point(18, 188);
            this.timeseries_form_button.Name = "timeseries_form_button";
            this.timeseries_form_button.Size = new System.Drawing.Size(142, 24);
            this.timeseries_form_button.TabIndex = 219;
            this.timeseries_form_button.Text = "timeseries outputs..";
            this.timeseries_form_button.UseVisualStyleBackColor = true;
            this.timeseries_form_button.Click += new System.EventHandler(this.timeseries_button_Click);
            // 
            // UTMgroupBox
            // 
            this.UTMgroupBox.Controls.Add(this.textBox6);
            this.UTMgroupBox.Controls.Add(this.UTMsouthcheck);
            this.UTMgroupBox.Controls.Add(this.UTMzonebox);
            this.UTMgroupBox.Location = new System.Drawing.Point(293, 131);
            this.UTMgroupBox.Name = "UTMgroupBox";
            this.UTMgroupBox.Size = new System.Drawing.Size(144, 65);
            this.UTMgroupBox.TabIndex = 218;
            this.UTMgroupBox.TabStop = false;
            this.UTMgroupBox.Visible = false;
            // 
            // textBox6
            // 
            this.textBox6.BorderStyle = System.Windows.Forms.BorderStyle.None;
            this.textBox6.Location = new System.Drawing.Point(6, 17);
            this.textBox6.Multiline = true;
            this.textBox6.Name = "textBox6";
            this.textBox6.ReadOnly = true;
            this.textBox6.Size = new System.Drawing.Size(87, 22);
            this.textBox6.TabIndex = 196;
            this.textBox6.Text = "UTM zone (1-60)";
            // 
            // UTMsouthcheck
            // 
            this.UTMsouthcheck.AutoSize = true;
            this.UTMsouthcheck.Location = new System.Drawing.Point(6, 42);
            this.UTMsouthcheck.Name = "UTMsouthcheck";
            this.UTMsouthcheck.Size = new System.Drawing.Size(128, 17);
            this.UTMsouthcheck.TabIndex = 197;
            this.UTMsouthcheck.Text = "Southern Hemisphere";
            this.UTMsouthcheck.UseVisualStyleBackColor = true;
            // 
            // UTMgridcheckbox
            // 
            this.UTMgridcheckbox.AutoSize = true;
            this.UTMgridcheckbox.Location = new System.Drawing.Point(293, 108);
            this.UTMgridcheckbox.Name = "UTMgridcheckbox";
            this.UTMgridcheckbox.Size = new System.Drawing.Size(82, 17);
            this.UTMgridcheckbox.TabIndex = 217;
            this.UTMgridcheckbox.Text = "Grid is UTM";
            this.UTMgridcheckbox.UseVisualStyleBackColor = true;
            this.UTMgridcheckbox.CheckedChanged += new System.EventHandler(this.UTMgridcheckbox_CheckedChanged);
            // 
            // textBox4
            // 
            this.textBox4.BorderStyle = System.Windows.Forms.BorderStyle.None;
            this.textBox4.Location = new System.Drawing.Point(293, 75);
            this.textBox4.Multiline = true;
            this.textBox4.Name = "textBox4";
            this.textBox4.ReadOnly = true;
            this.textBox4.Size = new System.Drawing.Size(138, 35);
            this.textBox4.TabIndex = 216;
            this.textBox4.Text = "only for British National Grid or UTM WGS84";
            // 
            // googleBeginDate
            // 
            this.googleBeginDate.AcceptsTab = true;
            this.googleBeginDate.Location = new System.Drawing.Point(165, 140);
            this.googleBeginDate.Name = "googleBeginDate";
            this.googleBeginDate.Size = new System.Drawing.Size(100, 20);
            this.googleBeginDate.TabIndex = 215;
            // 
            // googAnimationSaveInterval
            // 
            this.googAnimationSaveInterval.Location = new System.Drawing.Point(165, 114);
            this.googAnimationSaveInterval.Name = "googAnimationSaveInterval";
            this.googAnimationSaveInterval.Size = new System.Drawing.Size(56, 20);
            this.googAnimationSaveInterval.TabIndex = 212;
            this.googAnimationSaveInterval.Text = "1000";
            // 
            // saveintervalbox
            // 
            this.saveintervalbox.Location = new System.Drawing.Point(165, 56);
            this.saveintervalbox.Name = "saveintervalbox";
            this.saveintervalbox.Size = new System.Drawing.Size(56, 20);
            this.saveintervalbox.TabIndex = 201;
            this.saveintervalbox.Text = "1000";
            // 
            // label78
            // 
            this.label78.AutoSize = true;
            this.label78.Location = new System.Drawing.Point(29, 143);
            this.label78.Name = "label78";
            this.label78.Size = new System.Drawing.Size(120, 13);
            this.label78.TabIndex = 214;
            this.label78.Text = "begin date (yyyy-mm-dd)";
            // 
            // Run
            // 
            this.Run.Controls.Add(this.button4);
            this.Run.Controls.Add(this.version_lux_checkbox);
            this.Run.Controls.Add(this.groupBox2);
            this.Run.Controls.Add(this.calibration);
            this.Run.Controls.Add(this.Ik_ben_Marijn);
            this.Run.Controls.Add(this.groupBox7);
            this.Run.Location = new System.Drawing.Point(4, 22);
            this.Run.Name = "Run";
            this.Run.Size = new System.Drawing.Size(803, 293);
            this.Run.TabIndex = 8;
            this.Run.Text = "Run";
            this.Run.UseVisualStyleBackColor = true;
            // 
            // button4
            // 
            this.button4.Location = new System.Drawing.Point(102, 228);
            this.button4.Name = "button4";
            this.button4.Size = new System.Drawing.Size(163, 41);
            this.button4.TabIndex = 7;
            this.button4.Text = "now purely calculate terrain derivatives";
            this.button4.UseVisualStyleBackColor = true;
            this.button4.Click += new System.EventHandler(this.button4_Click);
            // 
            // version_lux_checkbox
            // 
            this.version_lux_checkbox.AutoSize = true;
            this.version_lux_checkbox.Location = new System.Drawing.Point(102, 167);
            this.version_lux_checkbox.Name = "version_lux_checkbox";
            this.version_lux_checkbox.Size = new System.Drawing.Size(115, 17);
            this.version_lux_checkbox.TabIndex = 6;
            this.version_lux_checkbox.Text = "Luxemburg version";
            this.version_lux_checkbox.UseVisualStyleBackColor = true;
            // 
            // groupBox2
            // 
            this.groupBox2.Controls.Add(this.label120);
            this.groupBox2.Controls.Add(this.calibration_ratio_reduction_parameter_textbox);
            this.groupBox2.Controls.Add(this.label119);
            this.groupBox2.Controls.Add(this.calibration_levels_textbox);
            this.groupBox2.Controls.Add(this.label116);
            this.groupBox2.Controls.Add(this.label118);
            this.groupBox2.Controls.Add(this.label117);
            this.groupBox2.Controls.Add(this.label115);
            this.groupBox2.Controls.Add(this.label114);
            this.groupBox2.Controls.Add(this.Sensitivity_button);
            this.groupBox2.Controls.Add(this.Calibration_button);
            this.groupBox2.Controls.Add(this.label113);
            this.groupBox2.Controls.Add(this.calibration_ratios_textbox);
            this.groupBox2.Location = new System.Drawing.Point(346, 32);
            this.groupBox2.Name = "groupBox2";
            this.groupBox2.Size = new System.Drawing.Size(422, 237);
            this.groupBox2.TabIndex = 4;
            this.groupBox2.TabStop = false;
            this.groupBox2.Text = "Calibration / sensitivity options";
            // 
            // label120
            // 
            this.label120.AutoSize = true;
            this.label120.Location = new System.Drawing.Point(39, 210);
            this.label120.Name = "label120";
            this.label120.Size = new System.Drawing.Size(199, 13);
            this.label120.TabIndex = 13;
            this.label120.Text = "1. describe the parameter values in code";
            // 
            // calibration_ratio_reduction_parameter_textbox
            // 
            this.calibration_ratio_reduction_parameter_textbox.Location = new System.Drawing.Point(338, 148);
            this.calibration_ratio_reduction_parameter_textbox.Name = "calibration_ratio_reduction_parameter_textbox";
            this.calibration_ratio_reduction_parameter_textbox.Size = new System.Drawing.Size(66, 20);
            this.calibration_ratio_reduction_parameter_textbox.TabIndex = 12;
            this.calibration_ratio_reduction_parameter_textbox.Text = "1.5";
            // 
            // label119
            // 
            this.label119.AutoSize = true;
            this.label119.Location = new System.Drawing.Point(39, 151);
            this.label119.Name = "label119";
            this.label119.Size = new System.Drawing.Size(208, 13);
            this.label119.TabIndex = 11;
            this.label119.Text = "5. reduction of variations per level (if smart)";
            // 
            // calibration_levels_textbox
            // 
            this.calibration_levels_textbox.Location = new System.Drawing.Point(338, 124);
            this.calibration_levels_textbox.Name = "calibration_levels_textbox";
            this.calibration_levels_textbox.Size = new System.Drawing.Size(66, 20);
            this.calibration_levels_textbox.TabIndex = 10;
            this.calibration_levels_textbox.Text = "3";
            // 
            // label116
            // 
            this.label116.AutoSize = true;
            this.label116.Location = new System.Drawing.Point(194, 35);
            this.label116.Name = "label116";
            this.label116.Size = new System.Drawing.Size(210, 13);
            this.label116.TabIndex = 9;
            this.label116.Text = "The optimal set of parameters will be stored";
            // 
            // label118
            // 
            this.label118.AutoSize = true;
            this.label118.Location = new System.Drawing.Point(39, 127);
            this.label118.Name = "label118";
            this.label118.Size = new System.Drawing.Size(97, 13);
            this.label118.TabIndex = 8;
            this.label118.Text = "4. levels (iterations)";
            // 
            // label117
            // 
            this.label117.AutoSize = true;
            this.label117.Location = new System.Drawing.Point(39, 80);
            this.label117.Name = "label117";
            this.label117.Size = new System.Drawing.Size(207, 13);
            this.label117.TabIndex = 7;
            this.label117.Text = "2. describe parameters to calibrate in code";
            // 
            // label115
            // 
            this.label115.AutoSize = true;
            this.label115.Location = new System.Drawing.Point(39, 56);
            this.label115.Name = "label115";
            this.label115.Size = new System.Drawing.Size(191, 13);
            this.label115.TabIndex = 5;
            this.label115.Text = "1. define the objective function in code";
            // 
            // label114
            // 
            this.label114.AutoSize = true;
            this.label114.Location = new System.Drawing.Point(39, 102);
            this.label114.Name = "label114";
            this.label114.Size = new System.Drawing.Size(132, 13);
            this.label114.TabIndex = 4;
            this.label114.Text = "3. variations per parameter";
            // 
            // Sensitivity_button
            // 
            this.Sensitivity_button.AutoSize = true;
            this.Sensitivity_button.Location = new System.Drawing.Point(22, 179);
            this.Sensitivity_button.Name = "Sensitivity_button";
            this.Sensitivity_button.Size = new System.Drawing.Size(200, 17);
            this.Sensitivity_button.TabIndex = 3;
            this.Sensitivity_button.Text = "Run sensitivity analysis (non-iterative)";
            this.Sensitivity_button.UseVisualStyleBackColor = true;
            this.Sensitivity_button.CheckedChanged += new System.EventHandler(this.radioButton2_CheckedChanged);
            // 
            // Calibration_button
            // 
            this.Calibration_button.AutoSize = true;
            this.Calibration_button.Location = new System.Drawing.Point(22, 33);
            this.Calibration_button.Name = "Calibration_button";
            this.Calibration_button.Size = new System.Drawing.Size(142, 17);
            this.Calibration_button.TabIndex = 2;
            this.Calibration_button.Text = "Run calibration (iterative)";
            this.Calibration_button.UseVisualStyleBackColor = true;
            this.Calibration_button.CheckedChanged += new System.EventHandler(this.radioButton1_CheckedChanged);
            // 
            // label113
            // 
            this.label113.AutoSize = true;
            this.label113.Location = new System.Drawing.Point(73, 39);
            this.label113.Name = "label113";
            this.label113.Size = new System.Drawing.Size(0, 13);
            this.label113.TabIndex = 1;
            // 
            // calibration_ratios_textbox
            // 
            this.calibration_ratios_textbox.Location = new System.Drawing.Point(218, 99);
            this.calibration_ratios_textbox.Name = "calibration_ratios_textbox";
            this.calibration_ratios_textbox.Size = new System.Drawing.Size(186, 20);
            this.calibration_ratios_textbox.TabIndex = 0;
            this.calibration_ratios_textbox.Text = "0.25;0.5;1;2;4";
            // 
            // calibration
            // 
            this.calibration.AutoSize = true;
            this.calibration.Location = new System.Drawing.Point(102, 190);
            this.calibration.Name = "calibration";
            this.calibration.Size = new System.Drawing.Size(125, 17);
            this.calibration.TabIndex = 5;
            this.calibration.Text = "Lessivage calibration";
            this.calibration.UseVisualStyleBackColor = true;
            this.calibration.CheckedChanged += new System.EventHandler(this.checkBox1_CheckedChanged_1);
            // 
            // Ik_ben_Marijn
            // 
            this.Ik_ben_Marijn.AutoSize = true;
            this.Ik_ben_Marijn.Location = new System.Drawing.Point(102, 146);
            this.Ik_ben_Marijn.Name = "Ik_ben_Marijn";
            this.Ik_ben_Marijn.Size = new System.Drawing.Size(87, 17);
            this.Ik_ben_Marijn.TabIndex = 4;
            this.Ik_ben_Marijn.Text = "Ik ben Marijn";
            this.Ik_ben_Marijn.UseVisualStyleBackColor = true;
            // 
            // groupBox7
            // 
            this.groupBox7.Controls.Add(this.runs_checkbox);
            this.groupBox7.Controls.Add(this.label16);
            this.groupBox7.Controls.Add(this.Number_runs_textbox);
            this.groupBox7.Location = new System.Drawing.Point(48, 32);
            this.groupBox7.Name = "groupBox7";
            this.groupBox7.Size = new System.Drawing.Size(277, 69);
            this.groupBox7.TabIndex = 3;
            this.groupBox7.TabStop = false;
            this.groupBox7.Text = "Please specify  the number of timesteps per run";
            // 
            // runs_checkbox
            // 
            this.runs_checkbox.AutoSize = true;
            this.runs_checkbox.Checked = true;
            this.runs_checkbox.Location = new System.Drawing.Point(54, 33);
            this.runs_checkbox.Name = "runs_checkbox";
            this.runs_checkbox.Size = new System.Drawing.Size(79, 17);
            this.runs_checkbox.TabIndex = 2;
            this.runs_checkbox.TabStop = true;
            this.runs_checkbox.Text = "runs (years)";
            this.runs_checkbox.UseVisualStyleBackColor = true;
            // 
            // label16
            // 
            this.label16.AutoSize = true;
            this.label16.Location = new System.Drawing.Point(73, 39);
            this.label16.Name = "label16";
            this.label16.Size = new System.Drawing.Size(0, 13);
            this.label16.TabIndex = 1;
            // 
            // Number_runs_textbox
            // 
            this.Number_runs_textbox.Location = new System.Drawing.Point(190, 30);
            this.Number_runs_textbox.Name = "Number_runs_textbox";
            this.Number_runs_textbox.Size = new System.Drawing.Size(55, 20);
            this.Number_runs_textbox.TabIndex = 0;
            this.Number_runs_textbox.Text = "1";
            // 
            // Input
            // 
            this.Input.Controls.Add(this.check_time_T);
            this.Input.Controls.Add(this.label99);
            this.Input.Controls.Add(this.label98);
            this.Input.Controls.Add(this.temp_input_filename_textbox);
            this.Input.Controls.Add(this.temp_constant_value_box);
            this.Input.Controls.Add(this.soil_specify_button);
            this.Input.Controls.Add(this.label88);
            this.Input.Controls.Add(this.groupBox3);
            this.Input.Controls.Add(this.explain_input_button);
            this.Input.Controls.Add(this.check_time_evap);
            this.Input.Controls.Add(this.check_time_infil);
            this.Input.Controls.Add(this.check_time_rain);
            this.Input.Controls.Add(this.check_time_till_fields);
            this.Input.Controls.Add(this.check_time_landuse);
            this.Input.Controls.Add(this.label29);
            this.Input.Controls.Add(this.check_space_DTM);
            this.Input.Controls.Add(this.groupBox13);
            this.Input.Controls.Add(this.tillfields_constant_textbox);
            this.Input.Controls.Add(this.tillfields_input_filename_textbox);
            this.Input.Controls.Add(this.evap_constant_value_box);
            this.Input.Controls.Add(this.evap_input_filename_textbox);
            this.Input.Controls.Add(this.infil_constant_value_box);
            this.Input.Controls.Add(this.infil_input_filename_textbox);
            this.Input.Controls.Add(this.rainfall_constant_value_box);
            this.Input.Controls.Add(this.landuse_constant_value_box);
            this.Input.Controls.Add(this.soildepth_constant_value_box);
            this.Input.Controls.Add(this.landuse_input_filename_textbox);
            this.Input.Controls.Add(this.soildepth_input_filename_textbox);
            this.Input.Controls.Add(this.rain_input_filename_textbox);
            this.Input.Controls.Add(this.dtm_input_filename_textbox);
            this.Input.Controls.Add(this.groupBox8);
            this.Input.Controls.Add(this.check_space_evap);
            this.Input.Controls.Add(this.check_space_infil);
            this.Input.Controls.Add(this.check_space_rain);
            this.Input.Controls.Add(this.check_space_till_fields);
            this.Input.Controls.Add(this.check_space_landuse);
            this.Input.Controls.Add(this.check_space_soildepth);
            this.Input.Controls.Add(this.label17);
            this.Input.Controls.Add(this.label15);
            this.Input.Controls.Add(this.label14);
            this.Input.Controls.Add(this.label7);
            this.Input.Controls.Add(label6);
            this.Input.Controls.Add(this.label5);
            this.Input.Controls.Add(this.label4);
            this.Input.Controls.Add(this.label3);
            this.Input.Controls.Add(this.label25);
            this.Input.Controls.Add(this.label23);
            this.Input.Location = new System.Drawing.Point(4, 22);
            this.Input.Name = "Input";
            this.Input.Size = new System.Drawing.Size(803, 293);
            this.Input.TabIndex = 0;
            this.Input.Text = "Inputs";
            this.Input.UseVisualStyleBackColor = true;
            // 
            // check_time_T
            // 
            this.check_time_T.AutoSize = true;
            this.check_time_T.Enabled = false;
            this.check_time_T.Location = new System.Drawing.Point(220, 251);
            this.check_time_T.Name = "check_time_T";
            this.check_time_T.Size = new System.Drawing.Size(15, 14);
            this.check_time_T.TabIndex = 154;
            this.check_time_T.UseVisualStyleBackColor = true;
            // 
            // temp_input_filename_textbox
            // 
            this.temp_input_filename_textbox.Enabled = false;
            this.temp_input_filename_textbox.Location = new System.Drawing.Point(258, 245);
            this.temp_input_filename_textbox.Name = "temp_input_filename_textbox";
            this.temp_input_filename_textbox.Size = new System.Drawing.Size(120, 20);
            this.temp_input_filename_textbox.TabIndex = 151;
            this.temp_input_filename_textbox.Text = "..";
            this.temp_input_filename_textbox.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // temp_constant_value_box
            // 
            this.temp_constant_value_box.Enabled = false;
            this.temp_constant_value_box.Location = new System.Drawing.Point(410, 245);
            this.temp_constant_value_box.Name = "temp_constant_value_box";
            this.temp_constant_value_box.Size = new System.Drawing.Size(120, 20);
            this.temp_constant_value_box.TabIndex = 150;
            this.temp_constant_value_box.Text = "10";
            this.temp_constant_value_box.TextChanged += new System.EventHandler(this.textBox3_TextChanged_1);
            // 
            // soil_specify_button
            // 
            this.soil_specify_button.Location = new System.Drawing.Point(188, 95);
            this.soil_specify_button.Name = "soil_specify_button";
            this.soil_specify_button.Size = new System.Drawing.Size(63, 20);
            this.soil_specify_button.TabIndex = 149;
            this.soil_specify_button.Text = "specify ..";
            this.soil_specify_button.UseVisualStyleBackColor = true;
            this.soil_specify_button.Click += new System.EventHandler(this.soil_specify_button_Click);
            // 
            // label88
            // 
            this.label88.Location = new System.Drawing.Point(19, 95);
            this.label88.Name = "label88";
            this.label88.Size = new System.Drawing.Size(104, 24);
            this.label88.TabIndex = 148;
            this.label88.Text = "soilproperties";
            this.label88.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            // 
            // explain_input_button
            // 
            this.explain_input_button.Location = new System.Drawing.Point(315, 25);
            this.explain_input_button.Name = "explain_input_button";
            this.explain_input_button.Size = new System.Drawing.Size(63, 20);
            this.explain_input_button.TabIndex = 147;
            this.explain_input_button.Text = "explain..";
            this.explain_input_button.UseVisualStyleBackColor = true;
            this.explain_input_button.Click += new System.EventHandler(this.button8_Click);
            // 
            // check_time_evap
            // 
            this.check_time_evap.AutoSize = true;
            this.check_time_evap.Location = new System.Drawing.Point(220, 222);
            this.check_time_evap.Name = "check_time_evap";
            this.check_time_evap.Size = new System.Drawing.Size(15, 14);
            this.check_time_evap.TabIndex = 145;
            this.check_time_evap.UseVisualStyleBackColor = true;
            this.check_time_evap.CheckedChanged += new System.EventHandler(this.check_time_evap_CheckedChanged);
            // 
            // check_time_infil
            // 
            this.check_time_infil.AutoSize = true;
            this.check_time_infil.Location = new System.Drawing.Point(220, 198);
            this.check_time_infil.Name = "check_time_infil";
            this.check_time_infil.Size = new System.Drawing.Size(15, 14);
            this.check_time_infil.TabIndex = 144;
            this.check_time_infil.UseVisualStyleBackColor = true;
            this.check_time_infil.CheckedChanged += new System.EventHandler(this.check_time_infil_CheckedChanged);
            // 
            // check_time_rain
            // 
            this.check_time_rain.AutoSize = true;
            this.check_time_rain.Location = new System.Drawing.Point(220, 174);
            this.check_time_rain.Name = "check_time_rain";
            this.check_time_rain.Size = new System.Drawing.Size(15, 14);
            this.check_time_rain.TabIndex = 143;
            this.check_time_rain.UseVisualStyleBackColor = true;
            this.check_time_rain.CheckedChanged += new System.EventHandler(this.check_time_rain_CheckedChanged);
            // 
            // check_time_till_fields
            // 
            this.check_time_till_fields.AutoSize = true;
            this.check_time_till_fields.Location = new System.Drawing.Point(220, 150);
            this.check_time_till_fields.Name = "check_time_till_fields";
            this.check_time_till_fields.Size = new System.Drawing.Size(15, 14);
            this.check_time_till_fields.TabIndex = 142;
            this.check_time_till_fields.UseVisualStyleBackColor = true;
            this.check_time_till_fields.CheckedChanged += new System.EventHandler(this.check_time_tillage_CheckedChanged);
            // 
            // check_time_landuse
            // 
            this.check_time_landuse.AutoSize = true;
            this.check_time_landuse.Location = new System.Drawing.Point(220, 123);
            this.check_time_landuse.Name = "check_time_landuse";
            this.check_time_landuse.Size = new System.Drawing.Size(15, 14);
            this.check_time_landuse.TabIndex = 141;
            this.check_time_landuse.UseVisualStyleBackColor = true;
            this.check_time_landuse.CheckedChanged += new System.EventHandler(this.check_time_landuse_CheckedChanged);
            // 
            // check_space_DTM
            // 
            this.check_space_DTM.AutoSize = true;
            this.check_space_DTM.Checked = true;
            this.check_space_DTM.CheckState = System.Windows.Forms.CheckState.Checked;
            this.check_space_DTM.Enabled = false;
            this.check_space_DTM.Location = new System.Drawing.Point(188, 51);
            this.check_space_DTM.Name = "check_space_DTM";
            this.check_space_DTM.Size = new System.Drawing.Size(15, 14);
            this.check_space_DTM.TabIndex = 138;
            this.check_space_DTM.UseVisualStyleBackColor = true;
            // 
            // tillfields_constant_textbox
            // 
            this.tillfields_constant_textbox.Location = new System.Drawing.Point(410, 145);
            this.tillfields_constant_textbox.Name = "tillfields_constant_textbox";
            this.tillfields_constant_textbox.ReadOnly = true;
            this.tillfields_constant_textbox.Size = new System.Drawing.Size(120, 20);
            this.tillfields_constant_textbox.TabIndex = 123;
            this.tillfields_constant_textbox.Text = "1";
            // 
            // tillfields_input_filename_textbox
            // 
            this.tillfields_input_filename_textbox.Enabled = false;
            this.tillfields_input_filename_textbox.Location = new System.Drawing.Point(258, 145);
            this.tillfields_input_filename_textbox.Name = "tillfields_input_filename_textbox";
            this.tillfields_input_filename_textbox.Size = new System.Drawing.Size(120, 20);
            this.tillfields_input_filename_textbox.TabIndex = 122;
            this.tillfields_input_filename_textbox.Text = "..";
            this.tillfields_input_filename_textbox.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            this.tillfields_input_filename_textbox.Click += new System.EventHandler(this.tillfields_input_filename_textbox_TextChanged);
            this.tillfields_input_filename_textbox.TextChanged += new System.EventHandler(this.tillfields_input_filename_textbox_TextChanged_1);
            // 
            // evap_constant_value_box
            // 
            this.evap_constant_value_box.Location = new System.Drawing.Point(410, 219);
            this.evap_constant_value_box.Name = "evap_constant_value_box";
            this.evap_constant_value_box.Size = new System.Drawing.Size(120, 20);
            this.evap_constant_value_box.TabIndex = 120;
            this.evap_constant_value_box.Text = "0.35";
            // 
            // evap_input_filename_textbox
            // 
            this.evap_input_filename_textbox.Enabled = false;
            this.evap_input_filename_textbox.Location = new System.Drawing.Point(258, 219);
            this.evap_input_filename_textbox.Name = "evap_input_filename_textbox";
            this.evap_input_filename_textbox.Size = new System.Drawing.Size(120, 20);
            this.evap_input_filename_textbox.TabIndex = 119;
            this.evap_input_filename_textbox.Text = "..";
            this.evap_input_filename_textbox.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            this.evap_input_filename_textbox.Click += new System.EventHandler(this.evap_input_filename_textbox_TextChanged);
            // 
            // infil_constant_value_box
            // 
            this.infil_constant_value_box.Location = new System.Drawing.Point(410, 193);
            this.infil_constant_value_box.Name = "infil_constant_value_box";
            this.infil_constant_value_box.Size = new System.Drawing.Size(120, 20);
            this.infil_constant_value_box.TabIndex = 117;
            this.infil_constant_value_box.Text = "0.150";
            // 
            // infil_input_filename_textbox
            // 
            this.infil_input_filename_textbox.Enabled = false;
            this.infil_input_filename_textbox.Location = new System.Drawing.Point(258, 193);
            this.infil_input_filename_textbox.Name = "infil_input_filename_textbox";
            this.infil_input_filename_textbox.Size = new System.Drawing.Size(120, 20);
            this.infil_input_filename_textbox.TabIndex = 116;
            this.infil_input_filename_textbox.Text = "..";
            this.infil_input_filename_textbox.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            this.infil_input_filename_textbox.Click += new System.EventHandler(this.infil_input_filename_textbox_TextChanged);
            // 
            // rainfall_constant_value_box
            // 
            this.rainfall_constant_value_box.Location = new System.Drawing.Point(410, 169);
            this.rainfall_constant_value_box.Name = "rainfall_constant_value_box";
            this.rainfall_constant_value_box.Size = new System.Drawing.Size(120, 20);
            this.rainfall_constant_value_box.TabIndex = 114;
            this.rainfall_constant_value_box.Text = "0.700";
            // 
            // landuse_constant_value_box
            // 
            this.landuse_constant_value_box.Location = new System.Drawing.Point(410, 120);
            this.landuse_constant_value_box.Name = "landuse_constant_value_box";
            this.landuse_constant_value_box.Size = new System.Drawing.Size(120, 20);
            this.landuse_constant_value_box.TabIndex = 113;
            this.landuse_constant_value_box.Text = "1";
            // 
            // soildepth_constant_value_box
            // 
            this.soildepth_constant_value_box.Location = new System.Drawing.Point(410, 76);
            this.soildepth_constant_value_box.Name = "soildepth_constant_value_box";
            this.soildepth_constant_value_box.Size = new System.Drawing.Size(120, 20);
            this.soildepth_constant_value_box.TabIndex = 112;
            this.soildepth_constant_value_box.Text = "100";
            // 
            // landuse_input_filename_textbox
            // 
            this.landuse_input_filename_textbox.Location = new System.Drawing.Point(258, 118);
            this.landuse_input_filename_textbox.Name = "landuse_input_filename_textbox";
            this.landuse_input_filename_textbox.Size = new System.Drawing.Size(120, 20);
            this.landuse_input_filename_textbox.TabIndex = 107;
            this.landuse_input_filename_textbox.Text = "..";
            this.landuse_input_filename_textbox.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            this.landuse_input_filename_textbox.Click += new System.EventHandler(this.landuse_input_filename_textbox_TextChanged);
            // 
            // soildepth_input_filename_textbox
            // 
            this.soildepth_input_filename_textbox.Location = new System.Drawing.Point(258, 76);
            this.soildepth_input_filename_textbox.Name = "soildepth_input_filename_textbox";
            this.soildepth_input_filename_textbox.Size = new System.Drawing.Size(120, 20);
            this.soildepth_input_filename_textbox.TabIndex = 105;
            this.soildepth_input_filename_textbox.Text = "..";
            this.soildepth_input_filename_textbox.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            this.soildepth_input_filename_textbox.Click += new System.EventHandler(this.soildepth_input_filename_textbox_TextChanged);
            // 
            // rain_input_filename_textbox
            // 
            this.rain_input_filename_textbox.Enabled = false;
            this.rain_input_filename_textbox.Location = new System.Drawing.Point(258, 169);
            this.rain_input_filename_textbox.Name = "rain_input_filename_textbox";
            this.rain_input_filename_textbox.Size = new System.Drawing.Size(120, 20);
            this.rain_input_filename_textbox.TabIndex = 103;
            this.rain_input_filename_textbox.Text = "..";
            this.rain_input_filename_textbox.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            this.rain_input_filename_textbox.Click += new System.EventHandler(this.rain_input_filename_textbox_TextChanged);
            this.rain_input_filename_textbox.TextChanged += new System.EventHandler(this.rain_input_filename_textbox_TextChanged_1);
            // 
            // dtm_input_filename_textbox
            // 
            this.dtm_input_filename_textbox.Location = new System.Drawing.Point(258, 50);
            this.dtm_input_filename_textbox.Name = "dtm_input_filename_textbox";
            this.dtm_input_filename_textbox.Size = new System.Drawing.Size(120, 20);
            this.dtm_input_filename_textbox.TabIndex = 100;
            this.dtm_input_filename_textbox.Text = "..";
            this.dtm_input_filename_textbox.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            this.dtm_input_filename_textbox.Click += new System.EventHandler(this.dtm_input_filename_textbox_Click);
            // 
            // groupBox8
            // 
            this.groupBox8.Controls.Add(this.fill_sinks_before_checkbox);
            this.groupBox8.Location = new System.Drawing.Point(589, 37);
            this.groupBox8.Name = "groupBox8";
            this.groupBox8.Size = new System.Drawing.Size(158, 48);
            this.groupBox8.TabIndex = 134;
            this.groupBox8.TabStop = false;
            this.groupBox8.Text = "before running: ";
            // 
            // check_space_evap
            // 
            this.check_space_evap.AutoSize = true;
            this.check_space_evap.Location = new System.Drawing.Point(188, 222);
            this.check_space_evap.Name = "check_space_evap";
            this.check_space_evap.Size = new System.Drawing.Size(15, 14);
            this.check_space_evap.TabIndex = 130;
            this.check_space_evap.UseVisualStyleBackColor = true;
            this.check_space_evap.CheckedChanged += new System.EventHandler(this.check_cnst_evap_CheckedChanged);
            // 
            // check_space_infil
            // 
            this.check_space_infil.AutoSize = true;
            this.check_space_infil.Location = new System.Drawing.Point(188, 198);
            this.check_space_infil.Name = "check_space_infil";
            this.check_space_infil.Size = new System.Drawing.Size(15, 14);
            this.check_space_infil.TabIndex = 129;
            this.check_space_infil.UseVisualStyleBackColor = true;
            this.check_space_infil.CheckedChanged += new System.EventHandler(this.check_cnst_infil_CheckedChanged);
            // 
            // check_space_rain
            // 
            this.check_space_rain.AutoSize = true;
            this.check_space_rain.Location = new System.Drawing.Point(188, 174);
            this.check_space_rain.Name = "check_space_rain";
            this.check_space_rain.Size = new System.Drawing.Size(15, 14);
            this.check_space_rain.TabIndex = 128;
            this.check_space_rain.UseVisualStyleBackColor = true;
            this.check_space_rain.CheckedChanged += new System.EventHandler(this.check_cnst_rain_CheckedChanged_1);
            // 
            // check_space_till_fields
            // 
            this.check_space_till_fields.AutoSize = true;
            this.check_space_till_fields.Location = new System.Drawing.Point(188, 150);
            this.check_space_till_fields.Name = "check_space_till_fields";
            this.check_space_till_fields.Size = new System.Drawing.Size(15, 14);
            this.check_space_till_fields.TabIndex = 127;
            this.check_space_till_fields.UseVisualStyleBackColor = true;
            this.check_space_till_fields.CheckedChanged += new System.EventHandler(this.check_cnst_till_fields_CheckedChanged);
            // 
            // check_space_landuse
            // 
            this.check_space_landuse.AutoSize = true;
            this.check_space_landuse.Location = new System.Drawing.Point(188, 123);
            this.check_space_landuse.Name = "check_space_landuse";
            this.check_space_landuse.Size = new System.Drawing.Size(15, 14);
            this.check_space_landuse.TabIndex = 126;
            this.check_space_landuse.UseVisualStyleBackColor = true;
            this.check_space_landuse.CheckedChanged += new System.EventHandler(this.check_cnst_landuse_CheckedChanged_1);
            // 
            // check_space_soildepth
            // 
            this.check_space_soildepth.AutoSize = true;
            this.check_space_soildepth.Location = new System.Drawing.Point(188, 76);
            this.check_space_soildepth.Name = "check_space_soildepth";
            this.check_space_soildepth.Size = new System.Drawing.Size(15, 14);
            this.check_space_soildepth.TabIndex = 125;
            this.check_space_soildepth.UseVisualStyleBackColor = true;
            this.check_space_soildepth.CheckedChanged += new System.EventHandler(this.check_cnst_soildepth_CheckedChanged_1);
            // 
            // label4
            // 
            this.label4.Location = new System.Drawing.Point(19, 116);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(104, 24);
            this.label4.TabIndex = 106;
            this.label4.Text = "landuse (classes)";
            this.label4.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            // 
            // label3
            // 
            this.label3.Location = new System.Drawing.Point(19, 71);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(104, 24);
            this.label3.TabIndex = 104;
            this.label3.Text = "soildepth [m]";
            this.label3.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            // 
            // label23
            // 
            this.label23.Location = new System.Drawing.Point(19, 45);
            this.label23.Name = "label23";
            this.label23.Size = new System.Drawing.Size(134, 24);
            this.label23.TabIndex = 56;
            this.label23.Text = "Digital Elevation Model [m]";
            this.label23.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            // 
            // Processes
            // 
            this.Processes.Controls.Add(this.Process_tabs);
            this.Processes.Location = new System.Drawing.Point(4, 22);
            this.Processes.Name = "Processes";
            this.Processes.Size = new System.Drawing.Size(803, 293);
            this.Processes.TabIndex = 6;
            this.Processes.Text = "Geomorphic processes";
            this.Processes.UseVisualStyleBackColor = true;
            // 
            // Process_tabs
            // 
            this.Process_tabs.Controls.Add(this.Water);
            this.Process_tabs.Controls.Add(this.Tillage);
            this.Process_tabs.Controls.Add(this.Creeper);
            this.Process_tabs.Controls.Add(Landsliding);
            this.Process_tabs.Controls.Add(this.Solifluction);
            this.Process_tabs.Controls.Add(this.Rock_weathering);
            this.Process_tabs.Controls.Add(this.Tectonics);
            this.Process_tabs.Controls.Add(this.treefall);
            this.Process_tabs.Location = new System.Drawing.Point(8, 14);
            this.Process_tabs.MaximumSize = new System.Drawing.Size(740, 276);
            this.Process_tabs.MinimumSize = new System.Drawing.Size(740, 276);
            this.Process_tabs.Name = "Process_tabs";
            this.Process_tabs.SelectedIndex = 0;
            this.Process_tabs.Size = new System.Drawing.Size(740, 276);
            this.Process_tabs.TabIndex = 0;
            // 
            // Water
            // 
            this.Water.Controls.Add(this.daily_water);
            this.Water.Controls.Add(this.label87);
            this.Water.Controls.Add(this.selectivity_constant_textbox);
            this.Water.Controls.Add(this.bio_protection_constant_textbox);
            this.Water.Controls.Add(this.erosion_threshold_textbox);
            this.Water.Controls.Add(this.rock_protection_constant_textbox);
            this.Water.Controls.Add(this.label90);
            this.Water.Controls.Add(this.label91);
            this.Water.Controls.Add(this.label92);
            this.Water.Controls.Add(this.parameter_n_textbox);
            this.Water.Controls.Add(this.parameter_conv_textbox);
            this.Water.Controls.Add(this.parameter_K_textbox);
            this.Water.Controls.Add(this.parameter_m_textbox);
            this.Water.Controls.Add(this.only_waterflow_checkbox);
            this.Water.Controls.Add(this.pictureBox1);
            this.Water.Controls.Add(this.label12);
            this.Water.Controls.Add(this.label11);
            this.Water.Controls.Add(this.label10);
            this.Water.Controls.Add(this.label9);
            this.Water.Controls.Add(this.Water_ero_checkbox);
            this.Water.Location = new System.Drawing.Point(4, 22);
            this.Water.Name = "Water";
            this.Water.Padding = new System.Windows.Forms.Padding(3);
            this.Water.Size = new System.Drawing.Size(732, 250);
            this.Water.TabIndex = 0;
            this.Water.Text = "Water erosion and deposition";
            this.Water.UseVisualStyleBackColor = true;
            // 
            // daily_water
            // 
            this.daily_water.AutoSize = true;
            this.daily_water.Location = new System.Drawing.Point(392, 16);
            this.daily_water.Name = "daily_water";
            this.daily_water.Size = new System.Drawing.Size(100, 17);
            this.daily_water.TabIndex = 29;
            this.daily_water.Text = "Daily water flow";
            this.daily_water.UseVisualStyleBackColor = true;
            this.daily_water.CheckedChanged += new System.EventHandler(this.daily_water_CheckedChanged);
            // 
            // label87
            // 
            this.label87.AutoSize = true;
            this.label87.Location = new System.Drawing.Point(101, 228);
            this.label87.Name = "label87";
            this.label87.Size = new System.Drawing.Size(136, 13);
            this.label87.TabIndex = 28;
            this.label87.Text = "selectivity change constant";
            // 
            // selectivity_constant_textbox
            // 
            this.selectivity_constant_textbox.Location = new System.Drawing.Point(26, 225);
            this.selectivity_constant_textbox.Name = "selectivity_constant_textbox";
            this.selectivity_constant_textbox.Size = new System.Drawing.Size(53, 20);
            this.selectivity_constant_textbox.TabIndex = 27;
            this.selectivity_constant_textbox.Text = "0";
            // 
            // bio_protection_constant_textbox
            // 
            this.bio_protection_constant_textbox.Location = new System.Drawing.Point(26, 199);
            this.bio_protection_constant_textbox.Name = "bio_protection_constant_textbox";
            this.bio_protection_constant_textbox.Size = new System.Drawing.Size(53, 20);
            this.bio_protection_constant_textbox.TabIndex = 21;
            this.bio_protection_constant_textbox.Text = "1";
            // 
            // erosion_threshold_textbox
            // 
            this.erosion_threshold_textbox.Location = new System.Drawing.Point(26, 147);
            this.erosion_threshold_textbox.Name = "erosion_threshold_textbox";
            this.erosion_threshold_textbox.Size = new System.Drawing.Size(53, 20);
            this.erosion_threshold_textbox.TabIndex = 20;
            this.erosion_threshold_textbox.Text = "0.01";
            // 
            // rock_protection_constant_textbox
            // 
            this.rock_protection_constant_textbox.Location = new System.Drawing.Point(26, 173);
            this.rock_protection_constant_textbox.Name = "rock_protection_constant_textbox";
            this.rock_protection_constant_textbox.Size = new System.Drawing.Size(53, 20);
            this.rock_protection_constant_textbox.TabIndex = 17;
            this.rock_protection_constant_textbox.Text = "1";
            // 
            // label90
            // 
            this.label90.AutoSize = true;
            this.label90.Location = new System.Drawing.Point(101, 150);
            this.label90.Name = "label90";
            this.label90.Size = new System.Drawing.Size(87, 13);
            this.label90.TabIndex = 24;
            this.label90.Text = "erosion threshold";
            // 
            // label91
            // 
            this.label91.AutoSize = true;
            this.label91.Location = new System.Drawing.Point(101, 202);
            this.label91.Name = "label91";
            this.label91.Size = new System.Drawing.Size(115, 13);
            this.label91.TabIndex = 23;
            this.label91.Text = "bio protection constant";
            // 
            // label92
            // 
            this.label92.AutoSize = true;
            this.label92.Location = new System.Drawing.Point(101, 176);
            this.label92.Name = "label92";
            this.label92.Size = new System.Drawing.Size(122, 13);
            this.label92.TabIndex = 22;
            this.label92.Text = "rock protection constant";
            // 
            // parameter_n_textbox
            // 
            this.parameter_n_textbox.Location = new System.Drawing.Point(26, 96);
            this.parameter_n_textbox.Name = "parameter_n_textbox";
            this.parameter_n_textbox.Size = new System.Drawing.Size(53, 20);
            this.parameter_n_textbox.TabIndex = 7;
            this.parameter_n_textbox.Text = "1.3";
            // 
            // parameter_conv_textbox
            // 
            this.parameter_conv_textbox.Location = new System.Drawing.Point(26, 47);
            this.parameter_conv_textbox.Name = "parameter_conv_textbox";
            this.parameter_conv_textbox.Size = new System.Drawing.Size(53, 20);
            this.parameter_conv_textbox.TabIndex = 6;
            this.parameter_conv_textbox.Text = "2";
            this.parameter_conv_textbox.TextChanged += new System.EventHandler(this.parameter_conv_textbox_TextChanged);
            // 
            // parameter_K_textbox
            // 
            this.parameter_K_textbox.Location = new System.Drawing.Point(26, 121);
            this.parameter_K_textbox.Name = "parameter_K_textbox";
            this.parameter_K_textbox.Size = new System.Drawing.Size(53, 20);
            this.parameter_K_textbox.TabIndex = 5;
            this.parameter_K_textbox.Text = "0.0003";
            // 
            // parameter_m_textbox
            // 
            this.parameter_m_textbox.Location = new System.Drawing.Point(26, 70);
            this.parameter_m_textbox.Name = "parameter_m_textbox";
            this.parameter_m_textbox.Size = new System.Drawing.Size(53, 20);
            this.parameter_m_textbox.TabIndex = 1;
            this.parameter_m_textbox.Text = "1.67";
            // 
            // only_waterflow_checkbox
            // 
            this.only_waterflow_checkbox.AutoSize = true;
            this.only_waterflow_checkbox.Location = new System.Drawing.Point(156, 16);
            this.only_waterflow_checkbox.Name = "only_waterflow_checkbox";
            this.only_waterflow_checkbox.Size = new System.Drawing.Size(219, 17);
            this.only_waterflow_checkbox.TabIndex = 14;
            this.only_waterflow_checkbox.Text = "Only calculate waterflow, no ero and dep";
            this.only_waterflow_checkbox.UseVisualStyleBackColor = true;
            // 
            // pictureBox1
            // 
            this.pictureBox1.Image = ((System.Drawing.Image)(resources.GetObject("pictureBox1.Image")));
            this.pictureBox1.Location = new System.Drawing.Point(266, 53);
            this.pictureBox1.Name = "pictureBox1";
            this.pictureBox1.Size = new System.Drawing.Size(180, 137);
            this.pictureBox1.TabIndex = 13;
            this.pictureBox1.TabStop = false;
            // 
            // label12
            // 
            this.label12.AutoSize = true;
            this.label12.Location = new System.Drawing.Point(101, 124);
            this.label12.Name = "label12";
            this.label12.Size = new System.Drawing.Size(66, 13);
            this.label12.TabIndex = 11;
            this.label12.Text = "K (erodibility)";
            // 
            // label11
            // 
            this.label11.AutoSize = true;
            this.label11.Location = new System.Drawing.Point(101, 50);
            this.label11.Name = "label11";
            this.label11.Size = new System.Drawing.Size(109, 13);
            this.label11.TabIndex = 10;
            this.label11.Text = "p (multiple flow factor)";
            this.label11.Click += new System.EventHandler(this.label11_Click);
            // 
            // label10
            // 
            this.label10.AutoSize = true;
            this.label10.Location = new System.Drawing.Point(101, 99);
            this.label10.Name = "label10";
            this.label10.Size = new System.Drawing.Size(106, 13);
            this.label10.TabIndex = 9;
            this.label10.Text = "n (exponent of slope)";
            this.label10.Click += new System.EventHandler(this.label10_Click);
            // 
            // label9
            // 
            this.label9.AutoSize = true;
            this.label9.Location = new System.Drawing.Point(101, 73);
            this.label9.Name = "label9";
            this.label9.Size = new System.Drawing.Size(146, 13);
            this.label9.TabIndex = 8;
            this.label9.Text = "m (exponent of overland flow)";
            // 
            // Water_ero_checkbox
            // 
            this.Water_ero_checkbox.AutoSize = true;
            this.Water_ero_checkbox.Checked = true;
            this.Water_ero_checkbox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.Water_ero_checkbox.Location = new System.Drawing.Point(26, 16);
            this.Water_ero_checkbox.Name = "Water_ero_checkbox";
            this.Water_ero_checkbox.Size = new System.Drawing.Size(124, 17);
            this.Water_ero_checkbox.TabIndex = 0;
            this.Water_ero_checkbox.Text = "Activate this process";
            this.Water_ero_checkbox.UseVisualStyleBackColor = true;
            this.Water_ero_checkbox.CheckedChanged += new System.EventHandler(this.Water_ero_checkbox_CheckedChanged);
            // 
            // Tillage
            // 
            this.Tillage.Controls.Add(this.pictureBox2);
            this.Tillage.Controls.Add(this.label20);
            this.Tillage.Controls.Add(this.trte);
            this.Tillage.Controls.Add(this.parameter_tillage_constant_textbox);
            this.Tillage.Controls.Add(this.parameter_ploughing_depth_textbox);
            this.Tillage.Controls.Add(this.Tillage_checkbox);
            this.Tillage.Location = new System.Drawing.Point(4, 22);
            this.Tillage.Name = "Tillage";
            this.Tillage.Padding = new System.Windows.Forms.Padding(3);
            this.Tillage.Size = new System.Drawing.Size(732, 250);
            this.Tillage.TabIndex = 1;
            this.Tillage.Text = "Tillage";
            this.Tillage.UseVisualStyleBackColor = true;
            // 
            // pictureBox2
            // 
            this.pictureBox2.Image = ((System.Drawing.Image)(resources.GetObject("pictureBox2.Image")));
            this.pictureBox2.Location = new System.Drawing.Point(276, 57);
            this.pictureBox2.Name = "pictureBox2";
            this.pictureBox2.Size = new System.Drawing.Size(180, 137);
            this.pictureBox2.TabIndex = 20;
            this.pictureBox2.TabStop = false;
            // 
            // label20
            // 
            this.label20.AutoSize = true;
            this.label20.Location = new System.Drawing.Point(128, 87);
            this.label20.Name = "label20";
            this.label20.Size = new System.Drawing.Size(78, 13);
            this.label20.TabIndex = 19;
            this.label20.Text = "tillage constant";
            // 
            // trte
            // 
            this.trte.AutoSize = true;
            this.trte.Location = new System.Drawing.Point(128, 61);
            this.trte.Name = "trte";
            this.trte.Size = new System.Drawing.Size(83, 13);
            this.trte.TabIndex = 18;
            this.trte.Text = "ploughing depth";
            // 
            // parameter_tillage_constant_textbox
            // 
            this.parameter_tillage_constant_textbox.Location = new System.Drawing.Point(53, 84);
            this.parameter_tillage_constant_textbox.Name = "parameter_tillage_constant_textbox";
            this.parameter_tillage_constant_textbox.Size = new System.Drawing.Size(53, 20);
            this.parameter_tillage_constant_textbox.TabIndex = 17;
            this.parameter_tillage_constant_textbox.Text = "0.08";
            // 
            // parameter_ploughing_depth_textbox
            // 
            this.parameter_ploughing_depth_textbox.AcceptsTab = true;
            this.parameter_ploughing_depth_textbox.Location = new System.Drawing.Point(53, 58);
            this.parameter_ploughing_depth_textbox.Name = "parameter_ploughing_depth_textbox";
            this.parameter_ploughing_depth_textbox.Size = new System.Drawing.Size(53, 20);
            this.parameter_ploughing_depth_textbox.TabIndex = 13;
            this.parameter_ploughing_depth_textbox.Text = "0.45";
            // 
            // Tillage_checkbox
            // 
            this.Tillage_checkbox.AutoSize = true;
            this.Tillage_checkbox.Location = new System.Drawing.Point(26, 16);
            this.Tillage_checkbox.Name = "Tillage_checkbox";
            this.Tillage_checkbox.Size = new System.Drawing.Size(124, 17);
            this.Tillage_checkbox.TabIndex = 1;
            this.Tillage_checkbox.Text = "Activate this process";
            this.Tillage_checkbox.UseVisualStyleBackColor = true;
            // 
            // Creeper
            // 
            this.Creeper.Controls.Add(this.creep_testing);
            this.Creeper.Controls.Add(this.pictureBox3);
            this.Creeper.Controls.Add(this.label19);
            this.Creeper.Controls.Add(this.parameter_diffusivity_textbox);
            this.Creeper.Controls.Add(this.creep_active_checkbox);
            this.Creeper.Location = new System.Drawing.Point(4, 22);
            this.Creeper.Name = "Creeper";
            this.Creeper.Size = new System.Drawing.Size(732, 250);
            this.Creeper.TabIndex = 6;
            this.Creeper.Text = "Creep";
            this.Creeper.UseVisualStyleBackColor = true;
            // 
            // creep_testing
            // 
            this.creep_testing.AutoSize = true;
            this.creep_testing.Location = new System.Drawing.Point(26, 108);
            this.creep_testing.Name = "creep_testing";
            this.creep_testing.Size = new System.Drawing.Size(88, 17);
            this.creep_testing.TabIndex = 26;
            this.creep_testing.Text = "Creep testing";
            this.creep_testing.UseVisualStyleBackColor = true;
            // 
            // pictureBox3
            // 
            this.pictureBox3.Image = ((System.Drawing.Image)(resources.GetObject("pictureBox3.Image")));
            this.pictureBox3.Location = new System.Drawing.Point(276, 57);
            this.pictureBox3.Name = "pictureBox3";
            this.pictureBox3.Size = new System.Drawing.Size(180, 137);
            this.pictureBox3.TabIndex = 25;
            this.pictureBox3.TabStop = false;
            // 
            // label19
            // 
            this.label19.AutoSize = true;
            this.label19.Location = new System.Drawing.Point(128, 63);
            this.label19.Name = "label19";
            this.label19.Size = new System.Drawing.Size(50, 13);
            this.label19.TabIndex = 23;
            this.label19.Text = "diffusivity";
            // 
            // parameter_diffusivity_textbox
            // 
            this.parameter_diffusivity_textbox.AcceptsTab = true;
            this.parameter_diffusivity_textbox.Location = new System.Drawing.Point(53, 60);
            this.parameter_diffusivity_textbox.Name = "parameter_diffusivity_textbox";
            this.parameter_diffusivity_textbox.Size = new System.Drawing.Size(53, 20);
            this.parameter_diffusivity_textbox.TabIndex = 21;
            this.parameter_diffusivity_textbox.Text = "0.05";
            // 
            // creep_active_checkbox
            // 
            this.creep_active_checkbox.AutoSize = true;
            this.creep_active_checkbox.Location = new System.Drawing.Point(26, 18);
            this.creep_active_checkbox.Name = "creep_active_checkbox";
            this.creep_active_checkbox.Size = new System.Drawing.Size(124, 17);
            this.creep_active_checkbox.TabIndex = 20;
            this.creep_active_checkbox.Text = "Activate this process";
            this.creep_active_checkbox.UseVisualStyleBackColor = true;
            // 
            // Solifluction
            // 
            this.Solifluction.Controls.Add(this.pictureBox5);
            this.Solifluction.Controls.Add(this.Solifluction_checkbox);
            this.Solifluction.Location = new System.Drawing.Point(4, 22);
            this.Solifluction.Name = "Solifluction";
            this.Solifluction.Size = new System.Drawing.Size(732, 250);
            this.Solifluction.TabIndex = 4;
            this.Solifluction.Text = "Solifluction";
            this.Solifluction.UseVisualStyleBackColor = true;
            // 
            // pictureBox5
            // 
            this.pictureBox5.Image = ((System.Drawing.Image)(resources.GetObject("pictureBox5.Image")));
            this.pictureBox5.Location = new System.Drawing.Point(276, 57);
            this.pictureBox5.Name = "pictureBox5";
            this.pictureBox5.Size = new System.Drawing.Size(180, 137);
            this.pictureBox5.TabIndex = 14;
            this.pictureBox5.TabStop = false;
            // 
            // Solifluction_checkbox
            // 
            this.Solifluction_checkbox.AutoSize = true;
            this.Solifluction_checkbox.Enabled = false;
            this.Solifluction_checkbox.Location = new System.Drawing.Point(36, 24);
            this.Solifluction_checkbox.Name = "Solifluction_checkbox";
            this.Solifluction_checkbox.Size = new System.Drawing.Size(124, 17);
            this.Solifluction_checkbox.TabIndex = 2;
            this.Solifluction_checkbox.Text = "Activate this process";
            this.Solifluction_checkbox.UseVisualStyleBackColor = true;
            // 
            // Rock_weathering
            // 
            this.Rock_weathering.Controls.Add(this.rockweath_method);
            this.Rock_weathering.Controls.Add(this.pictureBox6);
            this.Rock_weathering.Controls.Add(this.groupBox10);
            this.Rock_weathering.Controls.Add(this.groupBox9);
            this.Rock_weathering.Location = new System.Drawing.Point(4, 22);
            this.Rock_weathering.Name = "Rock_weathering";
            this.Rock_weathering.Size = new System.Drawing.Size(732, 250);
            this.Rock_weathering.TabIndex = 5;
            this.Rock_weathering.Text = "Rock weathering";
            this.Rock_weathering.UseVisualStyleBackColor = true;
            // 
            // rockweath_method
            // 
            this.rockweath_method.AllowDrop = true;
            this.rockweath_method.FormattingEnabled = true;
            this.rockweath_method.Items.AddRange(new object[] {
            "Humped",
            "Exponential (-P0 exp(-k1*dsoil))",
            "Function of infiltration (only with daily water flow)"});
            this.rockweath_method.Location = new System.Drawing.Point(26, 200);
            this.rockweath_method.Name = "rockweath_method";
            this.rockweath_method.Size = new System.Drawing.Size(121, 21);
            this.rockweath_method.TabIndex = 15;
            this.rockweath_method.Text = "Humped";
            this.rockweath_method.SelectedIndexChanged += new System.EventHandler(this.comboBox2_SelectedIndexChanged);
            // 
            // pictureBox6
            // 
            this.pictureBox6.Image = ((System.Drawing.Image)(resources.GetObject("pictureBox6.Image")));
            this.pictureBox6.Location = new System.Drawing.Point(276, 57);
            this.pictureBox6.Name = "pictureBox6";
            this.pictureBox6.Size = new System.Drawing.Size(180, 137);
            this.pictureBox6.TabIndex = 14;
            this.pictureBox6.TabStop = false;
            // 
            // groupBox10
            // 
            this.groupBox10.Controls.Add(this.Frost_weathering_checkbox);
            this.groupBox10.Enabled = false;
            this.groupBox10.Location = new System.Drawing.Point(250, 14);
            this.groupBox10.Name = "groupBox10";
            this.groupBox10.Size = new System.Drawing.Size(222, 179);
            this.groupBox10.TabIndex = 6;
            this.groupBox10.TabStop = false;
            this.groupBox10.Text = "Frost weathering ";
            this.groupBox10.Visible = false;
            // 
            // Frost_weathering_checkbox
            // 
            this.Frost_weathering_checkbox.AutoSize = true;
            this.Frost_weathering_checkbox.Enabled = false;
            this.Frost_weathering_checkbox.Location = new System.Drawing.Point(14, 19);
            this.Frost_weathering_checkbox.Name = "Frost_weathering_checkbox";
            this.Frost_weathering_checkbox.Size = new System.Drawing.Size(124, 17);
            this.Frost_weathering_checkbox.TabIndex = 3;
            this.Frost_weathering_checkbox.Text = "Activate this process";
            this.Frost_weathering_checkbox.UseVisualStyleBackColor = true;
            // 
            // Tectonics
            // 
            this.Tectonics.Controls.Add(this.groupBox14);
            this.Tectonics.Controls.Add(this.groupBox4);
            this.Tectonics.Location = new System.Drawing.Point(4, 22);
            this.Tectonics.Name = "Tectonics";
            this.Tectonics.Padding = new System.Windows.Forms.Padding(3);
            this.Tectonics.Size = new System.Drawing.Size(732, 250);
            this.Tectonics.TabIndex = 7;
            this.Tectonics.Text = "Tectonics";
            this.Tectonics.UseVisualStyleBackColor = true;
            // 
            // groupBox14
            // 
            this.groupBox14.Controls.Add(this.groupBox16);
            this.groupBox14.Controls.Add(this.Uplift_rate_textbox);
            this.groupBox14.Controls.Add(this.uplift_active_checkbox);
            this.groupBox14.Controls.Add(this.label39);
            this.groupBox14.Location = new System.Drawing.Point(176, 16);
            this.groupBox14.Name = "groupBox14";
            this.groupBox14.Size = new System.Drawing.Size(158, 209);
            this.groupBox14.TabIndex = 4;
            this.groupBox14.TabStop = false;
            this.groupBox14.Text = "Vertical uplift";
            // 
            // groupBox16
            // 
            this.groupBox16.Controls.Add(this.text_lift_col_less);
            this.groupBox16.Controls.Add(this.text_lift_col_more);
            this.groupBox16.Controls.Add(this.text_lift_row_less);
            this.groupBox16.Controls.Add(this.text_lift_row_more);
            this.groupBox16.Controls.Add(this.radio_lift_col_less_than);
            this.groupBox16.Controls.Add(this.radio_lift_row_more_than);
            this.groupBox16.Controls.Add(this.radio_lift_col_more_than);
            this.groupBox16.Controls.Add(this.radio_lift_row_less_than);
            this.groupBox16.Location = new System.Drawing.Point(13, 51);
            this.groupBox16.Name = "groupBox16";
            this.groupBox16.Size = new System.Drawing.Size(129, 105);
            this.groupBox16.TabIndex = 7;
            this.groupBox16.TabStop = false;
            this.groupBox16.Text = "For cells with:";
            // 
            // text_lift_col_less
            // 
            this.text_lift_col_less.Location = new System.Drawing.Point(63, 75);
            this.text_lift_col_less.Name = "text_lift_col_less";
            this.text_lift_col_less.Size = new System.Drawing.Size(54, 20);
            this.text_lift_col_less.TabIndex = 9;
            // 
            // text_lift_col_more
            // 
            this.text_lift_col_more.Location = new System.Drawing.Point(63, 56);
            this.text_lift_col_more.Name = "text_lift_col_more";
            this.text_lift_col_more.Size = new System.Drawing.Size(54, 20);
            this.text_lift_col_more.TabIndex = 8;
            // 
            // text_lift_row_less
            // 
            this.text_lift_row_less.Location = new System.Drawing.Point(63, 36);
            this.text_lift_row_less.Name = "text_lift_row_less";
            this.text_lift_row_less.Size = new System.Drawing.Size(54, 20);
            this.text_lift_row_less.TabIndex = 7;
            // 
            // text_lift_row_more
            // 
            this.text_lift_row_more.Location = new System.Drawing.Point(63, 16);
            this.text_lift_row_more.Name = "text_lift_row_more";
            this.text_lift_row_more.Size = new System.Drawing.Size(54, 20);
            this.text_lift_row_more.TabIndex = 6;
            // 
            // radio_lift_col_less_than
            // 
            this.radio_lift_col_less_than.AutoSize = true;
            this.radio_lift_col_less_than.Location = new System.Drawing.Point(6, 75);
            this.radio_lift_col_less_than.Name = "radio_lift_col_less_than";
            this.radio_lift_col_less_than.Size = new System.Drawing.Size(48, 17);
            this.radio_lift_col_less_than.TabIndex = 5;
            this.radio_lift_col_less_than.TabStop = true;
            this.radio_lift_col_less_than.Text = "col <";
            this.radio_lift_col_less_than.UseVisualStyleBackColor = true;
            // 
            // radio_lift_row_more_than
            // 
            this.radio_lift_row_more_than.AutoSize = true;
            this.radio_lift_row_more_than.Location = new System.Drawing.Point(6, 16);
            this.radio_lift_row_more_than.Name = "radio_lift_row_more_than";
            this.radio_lift_row_more_than.Size = new System.Drawing.Size(51, 17);
            this.radio_lift_row_more_than.TabIndex = 4;
            this.radio_lift_row_more_than.TabStop = true;
            this.radio_lift_row_more_than.Text = "row >";
            this.radio_lift_row_more_than.UseVisualStyleBackColor = true;
            // 
            // radio_lift_col_more_than
            // 
            this.radio_lift_col_more_than.AutoSize = true;
            this.radio_lift_col_more_than.Location = new System.Drawing.Point(6, 56);
            this.radio_lift_col_more_than.Name = "radio_lift_col_more_than";
            this.radio_lift_col_more_than.Size = new System.Drawing.Size(48, 17);
            this.radio_lift_col_more_than.TabIndex = 3;
            this.radio_lift_col_more_than.TabStop = true;
            this.radio_lift_col_more_than.Text = "col >";
            this.radio_lift_col_more_than.UseVisualStyleBackColor = true;
            // 
            // radio_lift_row_less_than
            // 
            this.radio_lift_row_less_than.AutoSize = true;
            this.radio_lift_row_less_than.Location = new System.Drawing.Point(6, 36);
            this.radio_lift_row_less_than.Name = "radio_lift_row_less_than";
            this.radio_lift_row_less_than.Size = new System.Drawing.Size(51, 17);
            this.radio_lift_row_less_than.TabIndex = 2;
            this.radio_lift_row_less_than.TabStop = true;
            this.radio_lift_row_less_than.Text = "row <";
            this.radio_lift_row_less_than.UseVisualStyleBackColor = true;
            // 
            // Uplift_rate_textbox
            // 
            this.Uplift_rate_textbox.Location = new System.Drawing.Point(13, 184);
            this.Uplift_rate_textbox.Name = "Uplift_rate_textbox";
            this.Uplift_rate_textbox.Size = new System.Drawing.Size(100, 20);
            this.Uplift_rate_textbox.TabIndex = 3;
            // 
            // uplift_active_checkbox
            // 
            this.uplift_active_checkbox.AutoSize = true;
            this.uplift_active_checkbox.Location = new System.Drawing.Point(13, 19);
            this.uplift_active_checkbox.Name = "uplift_active_checkbox";
            this.uplift_active_checkbox.Size = new System.Drawing.Size(65, 17);
            this.uplift_active_checkbox.TabIndex = 1;
            this.uplift_active_checkbox.Text = "Activate";
            this.uplift_active_checkbox.UseVisualStyleBackColor = true;
            // 
            // label39
            // 
            this.label39.AutoSize = true;
            this.label39.Location = new System.Drawing.Point(10, 168);
            this.label39.Name = "label39";
            this.label39.Size = new System.Drawing.Size(83, 13);
            this.label39.TabIndex = 2;
            this.label39.Text = "Uplift rate [m/a]:";
            // 
            // groupBox4
            // 
            this.groupBox4.Controls.Add(this.label38);
            this.groupBox4.Controls.Add(this.Tilting_rate_textbox);
            this.groupBox4.Controls.Add(this.groupBox15);
            this.groupBox4.Controls.Add(this.tilting_active_checkbox);
            this.groupBox4.Location = new System.Drawing.Point(13, 16);
            this.groupBox4.Name = "groupBox4";
            this.groupBox4.Size = new System.Drawing.Size(153, 210);
            this.groupBox4.TabIndex = 3;
            this.groupBox4.TabStop = false;
            this.groupBox4.Text = "Tilting";
            // 
            // label38
            // 
            this.label38.AutoSize = true;
            this.label38.Location = new System.Drawing.Point(9, 168);
            this.label38.Name = "label38";
            this.label38.Size = new System.Drawing.Size(114, 13);
            this.label38.TabIndex = 8;
            this.label38.Text = "Max alt change [m/a]: ";
            // 
            // Tilting_rate_textbox
            // 
            this.Tilting_rate_textbox.Location = new System.Drawing.Point(6, 184);
            this.Tilting_rate_textbox.Name = "Tilting_rate_textbox";
            this.Tilting_rate_textbox.Size = new System.Drawing.Size(100, 20);
            this.Tilting_rate_textbox.TabIndex = 7;
            // 
            // groupBox15
            // 
            this.groupBox15.Controls.Add(this.radio_tilt_col_max);
            this.groupBox15.Controls.Add(this.radio_tilt_row_zero);
            this.groupBox15.Controls.Add(this.radio_tilt_col_zero);
            this.groupBox15.Controls.Add(this.radio_tilt_row_max);
            this.groupBox15.Location = new System.Drawing.Point(6, 51);
            this.groupBox15.Name = "groupBox15";
            this.groupBox15.Size = new System.Drawing.Size(113, 105);
            this.groupBox15.TabIndex = 6;
            this.groupBox15.TabStop = false;
            this.groupBox15.Text = "Stability along:";
            // 
            // radio_tilt_col_max
            // 
            this.radio_tilt_col_max.AutoSize = true;
            this.radio_tilt_col_max.Location = new System.Drawing.Point(6, 79);
            this.radio_tilt_col_max.Name = "radio_tilt_col_max";
            this.radio_tilt_col_max.Size = new System.Drawing.Size(87, 17);
            this.radio_tilt_col_max.TabIndex = 5;
            this.radio_tilt_col_max.TabStop = true;
            this.radio_tilt_col_max.Text = "col = max col";
            this.radio_tilt_col_max.UseVisualStyleBackColor = true;
            // 
            // radio_tilt_row_zero
            // 
            this.radio_tilt_row_zero.AutoSize = true;
            this.radio_tilt_row_zero.Location = new System.Drawing.Point(6, 16);
            this.radio_tilt_row_zero.Name = "radio_tilt_row_zero";
            this.radio_tilt_row_zero.Size = new System.Drawing.Size(60, 17);
            this.radio_tilt_row_zero.TabIndex = 4;
            this.radio_tilt_row_zero.TabStop = true;
            this.radio_tilt_row_zero.Text = "row = 0";
            this.radio_tilt_row_zero.UseVisualStyleBackColor = true;
            // 
            // radio_tilt_col_zero
            // 
            this.radio_tilt_col_zero.AutoSize = true;
            this.radio_tilt_col_zero.Location = new System.Drawing.Point(6, 56);
            this.radio_tilt_col_zero.Name = "radio_tilt_col_zero";
            this.radio_tilt_col_zero.Size = new System.Drawing.Size(57, 17);
            this.radio_tilt_col_zero.TabIndex = 3;
            this.radio_tilt_col_zero.TabStop = true;
            this.radio_tilt_col_zero.Text = "col = 0";
            this.radio_tilt_col_zero.UseVisualStyleBackColor = true;
            // 
            // radio_tilt_row_max
            // 
            this.radio_tilt_row_max.AutoSize = true;
            this.radio_tilt_row_max.Location = new System.Drawing.Point(6, 36);
            this.radio_tilt_row_max.Name = "radio_tilt_row_max";
            this.radio_tilt_row_max.Size = new System.Drawing.Size(93, 17);
            this.radio_tilt_row_max.TabIndex = 2;
            this.radio_tilt_row_max.TabStop = true;
            this.radio_tilt_row_max.Text = "row = max row";
            this.radio_tilt_row_max.UseVisualStyleBackColor = true;
            // 
            // tilting_active_checkbox
            // 
            this.tilting_active_checkbox.AutoSize = true;
            this.tilting_active_checkbox.Location = new System.Drawing.Point(6, 19);
            this.tilting_active_checkbox.Name = "tilting_active_checkbox";
            this.tilting_active_checkbox.Size = new System.Drawing.Size(65, 17);
            this.tilting_active_checkbox.TabIndex = 0;
            this.tilting_active_checkbox.Text = "Activate";
            this.tilting_active_checkbox.UseVisualStyleBackColor = true;
            // 
            // treefall
            // 
            this.treefall.Controls.Add(this.tf_freq);
            this.treefall.Controls.Add(this.label112);
            this.treefall.Controls.Add(this.tf_age);
            this.treefall.Controls.Add(this.label111);
            this.treefall.Controls.Add(this.tf_growth);
            this.treefall.Controls.Add(this.label110);
            this.treefall.Controls.Add(this.tf_D);
            this.treefall.Controls.Add(this.label95);
            this.treefall.Controls.Add(this.label107);
            this.treefall.Controls.Add(this.tf_W);
            this.treefall.Controls.Add(this.treefall_checkbox);
            this.treefall.Location = new System.Drawing.Point(4, 22);
            this.treefall.Name = "treefall";
            this.treefall.Size = new System.Drawing.Size(732, 250);
            this.treefall.TabIndex = 8;
            this.treefall.Text = "Tree fall";
            this.treefall.UseVisualStyleBackColor = true;
            // 
            // tf_freq
            // 
            this.tf_freq.Location = new System.Drawing.Point(25, 162);
            this.tf_freq.Name = "tf_freq";
            this.tf_freq.Size = new System.Drawing.Size(53, 20);
            this.tf_freq.TabIndex = 30;
            this.tf_freq.Text = "0.00002";
            this.tf_freq.TextChanged += new System.EventHandler(this.textBox3_TextChanged_4);
            // 
            // label112
            // 
            this.label112.AutoSize = true;
            this.label112.Location = new System.Drawing.Point(100, 165);
            this.label112.Name = "label112";
            this.label112.Size = new System.Drawing.Size(132, 13);
            this.label112.TabIndex = 29;
            this.label112.Text = "fall frequency [trees/m2/a]";
            this.label112.Click += new System.EventHandler(this.label112_Click);
            // 
            // tf_age
            // 
            this.tf_age.Location = new System.Drawing.Point(25, 131);
            this.tf_age.Name = "tf_age";
            this.tf_age.Size = new System.Drawing.Size(53, 20);
            this.tf_age.TabIndex = 28;
            this.tf_age.Text = "300";
            // 
            // label111
            // 
            this.label111.AutoSize = true;
            this.label111.Location = new System.Drawing.Point(100, 134);
            this.label111.Name = "label111";
            this.label111.Size = new System.Drawing.Size(119, 13);
            this.label111.TabIndex = 27;
            this.label111.Text = "maximum age of tree [a]";
            this.label111.Click += new System.EventHandler(this.label111_Click);
            // 
            // tf_growth
            // 
            this.tf_growth.Location = new System.Drawing.Point(25, 103);
            this.tf_growth.Name = "tf_growth";
            this.tf_growth.Size = new System.Drawing.Size(53, 20);
            this.tf_growth.TabIndex = 26;
            this.tf_growth.Text = "150";
            // 
            // label110
            // 
            this.label110.AutoSize = true;
            this.label110.Location = new System.Drawing.Point(100, 106);
            this.label110.Name = "label110";
            this.label110.Size = new System.Drawing.Size(204, 13);
            this.label110.TabIndex = 25;
            this.label110.Text = "time it takes to reach these dimensions [a]";
            // 
            // tf_D
            // 
            this.tf_D.Location = new System.Drawing.Point(25, 77);
            this.tf_D.Name = "tf_D";
            this.tf_D.Size = new System.Drawing.Size(53, 20);
            this.tf_D.TabIndex = 24;
            this.tf_D.Text = "0.7";
            // 
            // label95
            // 
            this.label95.AutoSize = true;
            this.label95.Location = new System.Drawing.Point(100, 80);
            this.label95.Name = "label95";
            this.label95.Size = new System.Drawing.Size(145, 13);
            this.label95.TabIndex = 23;
            this.label95.Text = "maximum depth root mass [m]";
            this.label95.Click += new System.EventHandler(this.label95_Click_1);
            // 
            // label107
            // 
            this.label107.AutoSize = true;
            this.label107.Location = new System.Drawing.Point(100, 54);
            this.label107.Name = "label107";
            this.label107.Size = new System.Drawing.Size(158, 13);
            this.label107.TabIndex = 22;
            this.label107.Text = "maximum diameter root mass [m]";
            // 
            // tf_W
            // 
            this.tf_W.Location = new System.Drawing.Point(25, 51);
            this.tf_W.Name = "tf_W";
            this.tf_W.Size = new System.Drawing.Size(53, 20);
            this.tf_W.TabIndex = 21;
            this.tf_W.Text = "4";
            // 
            // treefall_checkbox
            // 
            this.treefall_checkbox.AutoSize = true;
            this.treefall_checkbox.Location = new System.Drawing.Point(25, 16);
            this.treefall_checkbox.Name = "treefall_checkbox";
            this.treefall_checkbox.Size = new System.Drawing.Size(124, 17);
            this.treefall_checkbox.TabIndex = 0;
            this.treefall_checkbox.Text = "Activate this process";
            this.treefall_checkbox.UseVisualStyleBackColor = true;
            this.treefall_checkbox.CheckedChanged += new System.EventHandler(this.checkBox1_CheckedChanged_2);
            // 
            // Creep_Checkbox
            // 
            this.Creep_Checkbox.AutoSize = true;
            this.Creep_Checkbox.Enabled = false;
            this.Creep_Checkbox.Location = new System.Drawing.Point(26, 16);
            this.Creep_Checkbox.Name = "Creep_Checkbox";
            this.Creep_Checkbox.Size = new System.Drawing.Size(124, 17);
            this.Creep_Checkbox.TabIndex = 1;
            this.Creep_Checkbox.Text = "Activate this process";
            this.Creep_Checkbox.UseVisualStyleBackColor = true;
            // 
            // tabControl1
            // 
            this.tabControl1.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.tabControl1.Controls.Add(this.Processes);
            this.tabControl1.Controls.Add(this.tabPage1);
            this.tabControl1.Controls.Add(this.tabPage2);
            this.tabControl1.Controls.Add(this.Input);
            this.tabControl1.Controls.Add(this.Run);
            this.tabControl1.Controls.Add(this.Output);
            this.tabControl1.Location = new System.Drawing.Point(4, 12);
            this.tabControl1.MaximumSize = new System.Drawing.Size(811, 319);
            this.tabControl1.MinimumSize = new System.Drawing.Size(811, 319);
            this.tabControl1.Name = "tabControl1";
            this.tabControl1.SelectedIndex = 0;
            this.tabControl1.Size = new System.Drawing.Size(811, 319);
            this.tabControl1.TabIndex = 143;
            // 
            // tabPage1
            // 
            this.tabPage1.Controls.Add(this.tabControl2);
            this.tabPage1.Location = new System.Drawing.Point(4, 22);
            this.tabPage1.Name = "tabPage1";
            this.tabPage1.Size = new System.Drawing.Size(803, 293);
            this.tabPage1.TabIndex = 9;
            this.tabPage1.Text = "Soil forming processes";
            this.tabPage1.UseVisualStyleBackColor = true;
            // 
            // tabControl2
            // 
            this.tabControl2.Controls.Add(this.physical);
            this.tabControl2.Controls.Add(this.chemical);
            this.tabControl2.Controls.Add(this.clay);
            this.tabControl2.Controls.Add(this.bioturbation);
            this.tabControl2.Controls.Add(this.carbon);
            this.tabControl2.Controls.Add(this.decalcification);
            this.tabControl2.Location = new System.Drawing.Point(16, 15);
            this.tabControl2.Name = "tabControl2";
            this.tabControl2.SelectedIndex = 0;
            this.tabControl2.Size = new System.Drawing.Size(759, 261);
            this.tabControl2.TabIndex = 0;
            // 
            // physical
            // 
            this.physical.Controls.Add(label49);
            this.physical.Controls.Add(label48);
            this.physical.Controls.Add(label47);
            this.physical.Controls.Add(label46);
            this.physical.Controls.Add(label45);
            this.physical.Controls.Add(label44);
            this.physical.Controls.Add(label43);
            this.physical.Controls.Add(label42);
            this.physical.Controls.Add(label41);
            this.physical.Controls.Add(this.upper_particle_fine_clay_textbox);
            this.physical.Controls.Add(this.upper_particle_clay_textbox);
            this.physical.Controls.Add(this.upper_particle_silt_textbox);
            this.physical.Controls.Add(this.upper_particle_sand_textbox);
            this.physical.Controls.Add(this.upper_particle_coarse_textbox);
            this.physical.Controls.Add(this.physical_weath_constant2);
            this.physical.Controls.Add(this.physical_weath_constant1);
            this.physical.Controls.Add(this.Physical_weath_C1_textbox);
            this.physical.Controls.Add(this.soil_phys_weath_checkbox);
            this.physical.Location = new System.Drawing.Point(4, 22);
            this.physical.Name = "physical";
            this.physical.Padding = new System.Windows.Forms.Padding(3);
            this.physical.Size = new System.Drawing.Size(751, 235);
            this.physical.TabIndex = 0;
            this.physical.Text = "Physical weathering";
            this.physical.UseVisualStyleBackColor = true;
            // 
            // upper_particle_fine_clay_textbox
            // 
            this.upper_particle_fine_clay_textbox.Location = new System.Drawing.Point(303, 147);
            this.upper_particle_fine_clay_textbox.Name = "upper_particle_fine_clay_textbox";
            this.upper_particle_fine_clay_textbox.Size = new System.Drawing.Size(100, 20);
            this.upper_particle_fine_clay_textbox.TabIndex = 9;
            this.upper_particle_fine_clay_textbox.Text = "0.0000001";
            // 
            // upper_particle_clay_textbox
            // 
            this.upper_particle_clay_textbox.Location = new System.Drawing.Point(303, 121);
            this.upper_particle_clay_textbox.Name = "upper_particle_clay_textbox";
            this.upper_particle_clay_textbox.Size = new System.Drawing.Size(100, 20);
            this.upper_particle_clay_textbox.TabIndex = 8;
            this.upper_particle_clay_textbox.Text = "0.000002";
            // 
            // upper_particle_silt_textbox
            // 
            this.upper_particle_silt_textbox.Location = new System.Drawing.Point(303, 95);
            this.upper_particle_silt_textbox.Name = "upper_particle_silt_textbox";
            this.upper_particle_silt_textbox.Size = new System.Drawing.Size(100, 20);
            this.upper_particle_silt_textbox.TabIndex = 7;
            this.upper_particle_silt_textbox.Text = "0.00005";
            // 
            // upper_particle_sand_textbox
            // 
            this.upper_particle_sand_textbox.Location = new System.Drawing.Point(303, 69);
            this.upper_particle_sand_textbox.Name = "upper_particle_sand_textbox";
            this.upper_particle_sand_textbox.Size = new System.Drawing.Size(100, 20);
            this.upper_particle_sand_textbox.TabIndex = 6;
            this.upper_particle_sand_textbox.Text = "0.002";
            // 
            // upper_particle_coarse_textbox
            // 
            this.upper_particle_coarse_textbox.Location = new System.Drawing.Point(303, 43);
            this.upper_particle_coarse_textbox.Name = "upper_particle_coarse_textbox";
            this.upper_particle_coarse_textbox.Size = new System.Drawing.Size(100, 20);
            this.upper_particle_coarse_textbox.TabIndex = 5;
            this.upper_particle_coarse_textbox.Text = "0.01";
            // 
            // physical_weath_constant2
            // 
            this.physical_weath_constant2.Location = new System.Drawing.Point(35, 95);
            this.physical_weath_constant2.Name = "physical_weath_constant2";
            this.physical_weath_constant2.Size = new System.Drawing.Size(100, 20);
            this.physical_weath_constant2.TabIndex = 4;
            this.physical_weath_constant2.Text = "5";
            // 
            // physical_weath_constant1
            // 
            this.physical_weath_constant1.Location = new System.Drawing.Point(35, 69);
            this.physical_weath_constant1.Name = "physical_weath_constant1";
            this.physical_weath_constant1.Size = new System.Drawing.Size(100, 20);
            this.physical_weath_constant1.TabIndex = 3;
            this.physical_weath_constant1.Text = "0.5";
            // 
            // Physical_weath_C1_textbox
            // 
            this.Physical_weath_C1_textbox.Location = new System.Drawing.Point(35, 43);
            this.Physical_weath_C1_textbox.Name = "Physical_weath_C1_textbox";
            this.Physical_weath_C1_textbox.Size = new System.Drawing.Size(100, 20);
            this.Physical_weath_C1_textbox.TabIndex = 2;
            this.Physical_weath_C1_textbox.Text = "0.000000004";
            // 
            // soil_phys_weath_checkbox
            // 
            this.soil_phys_weath_checkbox.AutoSize = true;
            this.soil_phys_weath_checkbox.Checked = true;
            this.soil_phys_weath_checkbox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.soil_phys_weath_checkbox.Location = new System.Drawing.Point(21, 6);
            this.soil_phys_weath_checkbox.Name = "soil_phys_weath_checkbox";
            this.soil_phys_weath_checkbox.Size = new System.Drawing.Size(124, 17);
            this.soil_phys_weath_checkbox.TabIndex = 1;
            this.soil_phys_weath_checkbox.Text = "Activate this process";
            this.soil_phys_weath_checkbox.UseVisualStyleBackColor = true;
            // 
            // chemical
            // 
            this.chemical.Controls.Add(label54);
            this.chemical.Controls.Add(label55);
            this.chemical.Controls.Add(label56);
            this.chemical.Controls.Add(label57);
            this.chemical.Controls.Add(label58);
            this.chemical.Controls.Add(label59);
            this.chemical.Controls.Add(this.specific_area_fine_clay_textbox);
            this.chemical.Controls.Add(this.specific_area_clay_textbox);
            this.chemical.Controls.Add(this.specific_area_silt_textbox);
            this.chemical.Controls.Add(this.specific_area_sand_textbox);
            this.chemical.Controls.Add(this.specific_area_coarse_textbox);
            this.chemical.Controls.Add(label53);
            this.chemical.Controls.Add(label50);
            this.chemical.Controls.Add(label51);
            this.chemical.Controls.Add(label52);
            this.chemical.Controls.Add(this.chem_weath_specific_coefficient_textbox);
            this.chemical.Controls.Add(this.chem_weath_depth_constant_textbox);
            this.chemical.Controls.Add(this.chem_weath_rate_constant_textbox);
            this.chemical.Controls.Add(this.soil_chem_weath_checkbox);
            this.chemical.Location = new System.Drawing.Point(4, 22);
            this.chemical.Name = "chemical";
            this.chemical.Padding = new System.Windows.Forms.Padding(3);
            this.chemical.Size = new System.Drawing.Size(751, 235);
            this.chemical.TabIndex = 1;
            this.chemical.Text = "Chemical weathering";
            this.chemical.UseVisualStyleBackColor = true;
            // 
            // specific_area_fine_clay_textbox
            // 
            this.specific_area_fine_clay_textbox.Location = new System.Drawing.Point(420, 142);
            this.specific_area_fine_clay_textbox.Name = "specific_area_fine_clay_textbox";
            this.specific_area_fine_clay_textbox.Size = new System.Drawing.Size(100, 20);
            this.specific_area_fine_clay_textbox.TabIndex = 24;
            this.specific_area_fine_clay_textbox.Text = "100000";
            // 
            // specific_area_clay_textbox
            // 
            this.specific_area_clay_textbox.Location = new System.Drawing.Point(420, 116);
            this.specific_area_clay_textbox.Name = "specific_area_clay_textbox";
            this.specific_area_clay_textbox.Size = new System.Drawing.Size(100, 20);
            this.specific_area_clay_textbox.TabIndex = 23;
            this.specific_area_clay_textbox.Text = "50000";
            // 
            // specific_area_silt_textbox
            // 
            this.specific_area_silt_textbox.Location = new System.Drawing.Point(420, 90);
            this.specific_area_silt_textbox.Name = "specific_area_silt_textbox";
            this.specific_area_silt_textbox.Size = new System.Drawing.Size(100, 20);
            this.specific_area_silt_textbox.TabIndex = 22;
            this.specific_area_silt_textbox.Text = "1000";
            // 
            // specific_area_sand_textbox
            // 
            this.specific_area_sand_textbox.Location = new System.Drawing.Point(420, 64);
            this.specific_area_sand_textbox.Name = "specific_area_sand_textbox";
            this.specific_area_sand_textbox.Size = new System.Drawing.Size(100, 20);
            this.specific_area_sand_textbox.TabIndex = 21;
            this.specific_area_sand_textbox.Text = "100";
            // 
            // specific_area_coarse_textbox
            // 
            this.specific_area_coarse_textbox.Location = new System.Drawing.Point(420, 38);
            this.specific_area_coarse_textbox.Name = "specific_area_coarse_textbox";
            this.specific_area_coarse_textbox.Size = new System.Drawing.Size(100, 20);
            this.specific_area_coarse_textbox.TabIndex = 20;
            this.specific_area_coarse_textbox.Text = "10";
            // 
            // chem_weath_specific_coefficient_textbox
            // 
            this.chem_weath_specific_coefficient_textbox.Location = new System.Drawing.Point(29, 90);
            this.chem_weath_specific_coefficient_textbox.Name = "chem_weath_specific_coefficient_textbox";
            this.chem_weath_specific_coefficient_textbox.Size = new System.Drawing.Size(100, 20);
            this.chem_weath_specific_coefficient_textbox.TabIndex = 15;
            this.chem_weath_specific_coefficient_textbox.Text = "1";
            // 
            // chem_weath_depth_constant_textbox
            // 
            this.chem_weath_depth_constant_textbox.Location = new System.Drawing.Point(29, 64);
            this.chem_weath_depth_constant_textbox.Name = "chem_weath_depth_constant_textbox";
            this.chem_weath_depth_constant_textbox.Size = new System.Drawing.Size(100, 20);
            this.chem_weath_depth_constant_textbox.TabIndex = 14;
            this.chem_weath_depth_constant_textbox.Text = "2.5";
            // 
            // chem_weath_rate_constant_textbox
            // 
            this.chem_weath_rate_constant_textbox.Location = new System.Drawing.Point(29, 38);
            this.chem_weath_rate_constant_textbox.Name = "chem_weath_rate_constant_textbox";
            this.chem_weath_rate_constant_textbox.Size = new System.Drawing.Size(100, 20);
            this.chem_weath_rate_constant_textbox.TabIndex = 13;
            this.chem_weath_rate_constant_textbox.Text = "0.000000004";
            // 
            // soil_chem_weath_checkbox
            // 
            this.soil_chem_weath_checkbox.AutoSize = true;
            this.soil_chem_weath_checkbox.Checked = true;
            this.soil_chem_weath_checkbox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.soil_chem_weath_checkbox.Location = new System.Drawing.Point(29, 6);
            this.soil_chem_weath_checkbox.Name = "soil_chem_weath_checkbox";
            this.soil_chem_weath_checkbox.Size = new System.Drawing.Size(124, 17);
            this.soil_chem_weath_checkbox.TabIndex = 1;
            this.soil_chem_weath_checkbox.Text = "Activate this process";
            this.soil_chem_weath_checkbox.UseVisualStyleBackColor = true;
            this.soil_chem_weath_checkbox.CheckedChanged += new System.EventHandler(this.soil_chem_weath_checkbox_CheckedChanged);
            // 
            // clay
            // 
            this.clay.Controls.Add(this.ct_Jagercikova);
            this.clay.Controls.Add(this.label109);
            this.clay.Controls.Add(this.label108);
            this.clay.Controls.Add(this.ct_dd_Jagercikova);
            this.clay.Controls.Add(this.ct_v0_Jagercikova);
            this.clay.Controls.Add(label13);
            this.clay.Controls.Add(this.ct_depth_decay);
            this.clay.Controls.Add(this.CT_depth_decay_checkbox);
            this.clay.Controls.Add(label69);
            this.clay.Controls.Add(label70);
            this.clay.Controls.Add(eluviation_rate_constant);
            this.clay.Controls.Add(this.eluviation_coefficient_textbox);
            this.clay.Controls.Add(this.maximum_eluviation_textbox);
            this.clay.Controls.Add(label72);
            this.clay.Controls.Add(label64);
            this.clay.Controls.Add(label65);
            this.clay.Controls.Add(label66);
            this.clay.Controls.Add(label67);
            this.clay.Controls.Add(this.clay_neoform_C2_textbox);
            this.clay.Controls.Add(this.clay_neoform_C1_textbox);
            this.clay.Controls.Add(this.clay_neoform_constant_textbox);
            this.clay.Controls.Add(label60);
            this.clay.Controls.Add(this.soil_clay_transloc_checkbox);
            this.clay.Location = new System.Drawing.Point(4, 22);
            this.clay.Name = "clay";
            this.clay.Size = new System.Drawing.Size(751, 235);
            this.clay.TabIndex = 2;
            this.clay.Text = "Clay dynamics";
            this.clay.UseVisualStyleBackColor = true;
            // 
            // ct_Jagercikova
            // 
            this.ct_Jagercikova.AutoSize = true;
            this.ct_Jagercikova.Location = new System.Drawing.Point(540, 52);
            this.ct_Jagercikova.Name = "ct_Jagercikova";
            this.ct_Jagercikova.Size = new System.Drawing.Size(179, 17);
            this.ct_Jagercikova.TabIndex = 62;
            this.ct_Jagercikova.Text = "Advection equation Jagercikova";
            this.ct_Jagercikova.UseVisualStyleBackColor = true;
            // 
            // label109
            // 
            this.label109.AutoSize = true;
            this.label109.Location = new System.Drawing.Point(602, 108);
            this.label109.Name = "label109";
            this.label109.Size = new System.Drawing.Size(98, 13);
            this.label109.TabIndex = 61;
            this.label109.Text = "depth decay [cm-1]";
            // 
            // label108
            // 
            this.label108.AutoSize = true;
            this.label108.Location = new System.Drawing.Point(597, 78);
            this.label108.Name = "label108";
            this.label108.Size = new System.Drawing.Size(148, 13);
            this.label108.TabIndex = 60;
            this.label108.Text = "surface advection v0 [cm a-1]";
            // 
            // ct_dd_Jagercikova
            // 
            this.ct_dd_Jagercikova.Location = new System.Drawing.Point(540, 104);
            this.ct_dd_Jagercikova.Name = "ct_dd_Jagercikova";
            this.ct_dd_Jagercikova.Size = new System.Drawing.Size(51, 20);
            this.ct_dd_Jagercikova.TabIndex = 58;
            this.ct_dd_Jagercikova.Text = "0.09";
            // 
            // ct_v0_Jagercikova
            // 
            this.ct_v0_Jagercikova.Location = new System.Drawing.Point(540, 75);
            this.ct_v0_Jagercikova.Name = "ct_v0_Jagercikova";
            this.ct_v0_Jagercikova.Size = new System.Drawing.Size(51, 20);
            this.ct_v0_Jagercikova.TabIndex = 57;
            this.ct_v0_Jagercikova.Text = "0.18";
            // 
            // ct_depth_decay
            // 
            this.ct_depth_decay.Location = new System.Drawing.Point(303, 169);
            this.ct_depth_decay.Name = "ct_depth_decay";
            this.ct_depth_decay.Size = new System.Drawing.Size(100, 20);
            this.ct_depth_decay.TabIndex = 55;
            this.ct_depth_decay.Text = "2";
            // 
            // CT_depth_decay_checkbox
            // 
            this.CT_depth_decay_checkbox.AutoSize = true;
            this.CT_depth_decay_checkbox.Checked = true;
            this.CT_depth_decay_checkbox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.CT_depth_decay_checkbox.Location = new System.Drawing.Point(304, 146);
            this.CT_depth_decay_checkbox.Name = "CT_depth_decay_checkbox";
            this.CT_depth_decay_checkbox.Size = new System.Drawing.Size(137, 17);
            this.CT_depth_decay_checkbox.TabIndex = 54;
            this.CT_depth_decay_checkbox.Text = "Depth decay constant?";
            this.CT_depth_decay_checkbox.UseVisualStyleBackColor = true;
            // 
            // eluviation_coefficient_textbox
            // 
            this.eluviation_coefficient_textbox.Location = new System.Drawing.Point(304, 101);
            this.eluviation_coefficient_textbox.Name = "eluviation_coefficient_textbox";
            this.eluviation_coefficient_textbox.Size = new System.Drawing.Size(100, 20);
            this.eluviation_coefficient_textbox.TabIndex = 49;
            this.eluviation_coefficient_textbox.Text = "2";
            // 
            // maximum_eluviation_textbox
            // 
            this.maximum_eluviation_textbox.Location = new System.Drawing.Point(304, 75);
            this.maximum_eluviation_textbox.Name = "maximum_eluviation_textbox";
            this.maximum_eluviation_textbox.Size = new System.Drawing.Size(100, 20);
            this.maximum_eluviation_textbox.TabIndex = 48;
            this.maximum_eluviation_textbox.Text = "0.007";
            // 
            // clay_neoform_C2_textbox
            // 
            this.clay_neoform_C2_textbox.Location = new System.Drawing.Point(25, 127);
            this.clay_neoform_C2_textbox.Name = "clay_neoform_C2_textbox";
            this.clay_neoform_C2_textbox.Size = new System.Drawing.Size(100, 20);
            this.clay_neoform_C2_textbox.TabIndex = 42;
            this.clay_neoform_C2_textbox.Text = "20";
            // 
            // clay_neoform_C1_textbox
            // 
            this.clay_neoform_C1_textbox.Location = new System.Drawing.Point(25, 101);
            this.clay_neoform_C1_textbox.Name = "clay_neoform_C1_textbox";
            this.clay_neoform_C1_textbox.Size = new System.Drawing.Size(100, 20);
            this.clay_neoform_C1_textbox.TabIndex = 41;
            this.clay_neoform_C1_textbox.Text = "1";
            // 
            // clay_neoform_constant_textbox
            // 
            this.clay_neoform_constant_textbox.Location = new System.Drawing.Point(25, 75);
            this.clay_neoform_constant_textbox.Name = "clay_neoform_constant_textbox";
            this.clay_neoform_constant_textbox.Size = new System.Drawing.Size(100, 20);
            this.clay_neoform_constant_textbox.TabIndex = 40;
            this.clay_neoform_constant_textbox.Text = "0.5";
            // 
            // soil_clay_transloc_checkbox
            // 
            this.soil_clay_transloc_checkbox.AutoSize = true;
            this.soil_clay_transloc_checkbox.Checked = true;
            this.soil_clay_transloc_checkbox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.soil_clay_transloc_checkbox.Location = new System.Drawing.Point(26, 12);
            this.soil_clay_transloc_checkbox.Name = "soil_clay_transloc_checkbox";
            this.soil_clay_transloc_checkbox.Size = new System.Drawing.Size(124, 17);
            this.soil_clay_transloc_checkbox.TabIndex = 1;
            this.soil_clay_transloc_checkbox.Text = "Activate this process";
            this.soil_clay_transloc_checkbox.UseVisualStyleBackColor = true;
            // 
            // bioturbation
            // 
            this.bioturbation.Controls.Add(label68);
            this.bioturbation.Controls.Add(label71);
            this.bioturbation.Controls.Add(label73);
            this.bioturbation.Controls.Add(this.bioturbation_depth_decay_textbox);
            this.bioturbation.Controls.Add(this.potential_bioturbation_textbox);
            this.bioturbation.Controls.Add(this.soil_bioturb_checkbox);
            this.bioturbation.Location = new System.Drawing.Point(4, 22);
            this.bioturbation.Name = "bioturbation";
            this.bioturbation.Size = new System.Drawing.Size(751, 235);
            this.bioturbation.TabIndex = 3;
            this.bioturbation.Text = "Bioturbation";
            this.bioturbation.UseVisualStyleBackColor = true;
            // 
            // bioturbation_depth_decay_textbox
            // 
            this.bioturbation_depth_decay_textbox.Location = new System.Drawing.Point(26, 74);
            this.bioturbation_depth_decay_textbox.Name = "bioturbation_depth_decay_textbox";
            this.bioturbation_depth_decay_textbox.Size = new System.Drawing.Size(100, 20);
            this.bioturbation_depth_decay_textbox.TabIndex = 56;
            this.bioturbation_depth_decay_textbox.Text = "2.5";
            // 
            // potential_bioturbation_textbox
            // 
            this.potential_bioturbation_textbox.Location = new System.Drawing.Point(26, 48);
            this.potential_bioturbation_textbox.Name = "potential_bioturbation_textbox";
            this.potential_bioturbation_textbox.Size = new System.Drawing.Size(100, 20);
            this.potential_bioturbation_textbox.TabIndex = 55;
            this.potential_bioturbation_textbox.Text = "6";
            // 
            // soil_bioturb_checkbox
            // 
            this.soil_bioturb_checkbox.AutoSize = true;
            this.soil_bioturb_checkbox.Checked = true;
            this.soil_bioturb_checkbox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.soil_bioturb_checkbox.Location = new System.Drawing.Point(26, 12);
            this.soil_bioturb_checkbox.Name = "soil_bioturb_checkbox";
            this.soil_bioturb_checkbox.Size = new System.Drawing.Size(124, 17);
            this.soil_bioturb_checkbox.TabIndex = 1;
            this.soil_bioturb_checkbox.Text = "Activate this process";
            this.soil_bioturb_checkbox.UseVisualStyleBackColor = true;
            // 
            // carbon
            // 
            this.carbon.Controls.Add(this.carbon_o_decomp_rate_textbox);
            this.carbon.Controls.Add(label86);
            this.carbon.Controls.Add(this.carbon_y_decomp_rate_textbox);
            this.carbon.Controls.Add(this.carbon_o_twi_decay_textbox);
            this.carbon.Controls.Add(label85);
            this.carbon.Controls.Add(this.carbon_y_twi_decay_textbox);
            this.carbon.Controls.Add(this.carbon_o_depth_decay_textbox);
            this.carbon.Controls.Add(label84);
            this.carbon.Controls.Add(this.carbon_y_depth_decay_textbox);
            this.carbon.Controls.Add(label83);
            this.carbon.Controls.Add(label82);
            this.carbon.Controls.Add(label80);
            this.carbon.Controls.Add(label77);
            this.carbon.Controls.Add(label81);
            this.carbon.Controls.Add(this.carbon_humification_fraction_textbox);
            this.carbon.Controls.Add(label74);
            this.carbon.Controls.Add(label75);
            this.carbon.Controls.Add(label76);
            this.carbon.Controls.Add(this.carbon_depth_decay_textbox);
            this.carbon.Controls.Add(this.carbon_input_textbox);
            this.carbon.Controls.Add(this.soil_carbon_cycle_checkbox);
            this.carbon.Location = new System.Drawing.Point(4, 22);
            this.carbon.Name = "carbon";
            this.carbon.Size = new System.Drawing.Size(751, 235);
            this.carbon.TabIndex = 4;
            this.carbon.Text = "Carbon Cycle";
            this.carbon.UseVisualStyleBackColor = true;
            // 
            // carbon_o_decomp_rate_textbox
            // 
            this.carbon_o_decomp_rate_textbox.Location = new System.Drawing.Point(453, 111);
            this.carbon_o_decomp_rate_textbox.Name = "carbon_o_decomp_rate_textbox";
            this.carbon_o_decomp_rate_textbox.Size = new System.Drawing.Size(100, 20);
            this.carbon_o_decomp_rate_textbox.TabIndex = 81;
            this.carbon_o_decomp_rate_textbox.Text = "0.005";
            // 
            // carbon_y_decomp_rate_textbox
            // 
            this.carbon_y_decomp_rate_textbox.Location = new System.Drawing.Point(347, 111);
            this.carbon_y_decomp_rate_textbox.Name = "carbon_y_decomp_rate_textbox";
            this.carbon_y_decomp_rate_textbox.Size = new System.Drawing.Size(100, 20);
            this.carbon_y_decomp_rate_textbox.TabIndex = 79;
            this.carbon_y_decomp_rate_textbox.Text = "0.01";
            // 
            // carbon_o_twi_decay_textbox
            // 
            this.carbon_o_twi_decay_textbox.Location = new System.Drawing.Point(453, 163);
            this.carbon_o_twi_decay_textbox.Name = "carbon_o_twi_decay_textbox";
            this.carbon_o_twi_decay_textbox.Size = new System.Drawing.Size(100, 20);
            this.carbon_o_twi_decay_textbox.TabIndex = 78;
            this.carbon_o_twi_decay_textbox.Text = "0.03";
            // 
            // carbon_y_twi_decay_textbox
            // 
            this.carbon_y_twi_decay_textbox.Location = new System.Drawing.Point(347, 163);
            this.carbon_y_twi_decay_textbox.Name = "carbon_y_twi_decay_textbox";
            this.carbon_y_twi_decay_textbox.Size = new System.Drawing.Size(100, 20);
            this.carbon_y_twi_decay_textbox.TabIndex = 76;
            this.carbon_y_twi_decay_textbox.Text = "0.03";
            // 
            // carbon_o_depth_decay_textbox
            // 
            this.carbon_o_depth_decay_textbox.Location = new System.Drawing.Point(453, 137);
            this.carbon_o_depth_decay_textbox.Name = "carbon_o_depth_decay_textbox";
            this.carbon_o_depth_decay_textbox.Size = new System.Drawing.Size(100, 20);
            this.carbon_o_depth_decay_textbox.TabIndex = 75;
            this.carbon_o_depth_decay_textbox.Text = "8";
            // 
            // carbon_y_depth_decay_textbox
            // 
            this.carbon_y_depth_decay_textbox.Location = new System.Drawing.Point(347, 137);
            this.carbon_y_depth_decay_textbox.Name = "carbon_y_depth_decay_textbox";
            this.carbon_y_depth_decay_textbox.Size = new System.Drawing.Size(100, 20);
            this.carbon_y_depth_decay_textbox.TabIndex = 73;
            this.carbon_y_depth_decay_textbox.Text = "8";
            // 
            // carbon_humification_fraction_textbox
            // 
            this.carbon_humification_fraction_textbox.Location = new System.Drawing.Point(23, 117);
            this.carbon_humification_fraction_textbox.Name = "carbon_humification_fraction_textbox";
            this.carbon_humification_fraction_textbox.Size = new System.Drawing.Size(100, 20);
            this.carbon_humification_fraction_textbox.TabIndex = 65;
            this.carbon_humification_fraction_textbox.Text = "0.8";
            // 
            // carbon_depth_decay_textbox
            // 
            this.carbon_depth_decay_textbox.Location = new System.Drawing.Point(23, 88);
            this.carbon_depth_decay_textbox.Name = "carbon_depth_decay_textbox";
            this.carbon_depth_decay_textbox.Size = new System.Drawing.Size(100, 20);
            this.carbon_depth_decay_textbox.TabIndex = 61;
            this.carbon_depth_decay_textbox.Text = "8";
            // 
            // carbon_input_textbox
            // 
            this.carbon_input_textbox.Location = new System.Drawing.Point(23, 62);
            this.carbon_input_textbox.Name = "carbon_input_textbox";
            this.carbon_input_textbox.Size = new System.Drawing.Size(100, 20);
            this.carbon_input_textbox.TabIndex = 60;
            this.carbon_input_textbox.Text = "1.5";
            // 
            // soil_carbon_cycle_checkbox
            // 
            this.soil_carbon_cycle_checkbox.AutoSize = true;
            this.soil_carbon_cycle_checkbox.Checked = true;
            this.soil_carbon_cycle_checkbox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.soil_carbon_cycle_checkbox.Location = new System.Drawing.Point(25, 14);
            this.soil_carbon_cycle_checkbox.Name = "soil_carbon_cycle_checkbox";
            this.soil_carbon_cycle_checkbox.Size = new System.Drawing.Size(124, 17);
            this.soil_carbon_cycle_checkbox.TabIndex = 2;
            this.soil_carbon_cycle_checkbox.Text = "Activate this process";
            this.soil_carbon_cycle_checkbox.UseVisualStyleBackColor = true;
            // 
            // decalcification
            // 
            this.decalcification.Controls.Add(this.label94);
            this.decalcification.Controls.Add(this.ini_CaCO3_content);
            this.decalcification.Controls.Add(this.decalcification_checkbox);
            this.decalcification.Location = new System.Drawing.Point(4, 22);
            this.decalcification.Name = "decalcification";
            this.decalcification.Size = new System.Drawing.Size(751, 235);
            this.decalcification.TabIndex = 5;
            this.decalcification.Text = "Decalcification";
            this.decalcification.UseVisualStyleBackColor = true;
            // 
            // label94
            // 
            this.label94.AutoSize = true;
            this.label94.Location = new System.Drawing.Point(138, 42);
            this.label94.Name = "label94";
            this.label94.Size = new System.Drawing.Size(107, 13);
            this.label94.TabIndex = 2;
            this.label94.Text = "Initial CaCO3 content";
            // 
            // ini_CaCO3_content
            // 
            this.ini_CaCO3_content.Location = new System.Drawing.Point(32, 39);
            this.ini_CaCO3_content.Name = "ini_CaCO3_content";
            this.ini_CaCO3_content.Size = new System.Drawing.Size(100, 20);
            this.ini_CaCO3_content.TabIndex = 1;
            this.ini_CaCO3_content.Text = "0.1";
            // 
            // decalcification_checkbox
            // 
            this.decalcification_checkbox.AutoSize = true;
            this.decalcification_checkbox.Location = new System.Drawing.Point(32, 15);
            this.decalcification_checkbox.Name = "decalcification_checkbox";
            this.decalcification_checkbox.Size = new System.Drawing.Size(124, 17);
            this.decalcification_checkbox.TabIndex = 0;
            this.decalcification_checkbox.Text = "Activate this process";
            this.decalcification_checkbox.UseVisualStyleBackColor = true;
            this.decalcification_checkbox.CheckedChanged += new System.EventHandler(this.decalcification_checkbox_CheckedChanged);
            // 
            // tabPage2
            // 
            this.tabPage2.Controls.Add(this.check_scaling_daily_weather);
            this.tabPage2.Controls.Add(this.label106);
            this.tabPage2.Controls.Add(this.snow_threshold_textbox);
            this.tabPage2.Controls.Add(this.label105);
            this.tabPage2.Controls.Add(this.snowmelt_factor_textbox);
            this.tabPage2.Controls.Add(this.label104);
            this.tabPage2.Controls.Add(this.latitude_min);
            this.tabPage2.Controls.Add(this.label103);
            this.tabPage2.Controls.Add(this.latitude_deg);
            this.tabPage2.Controls.Add(this.label100);
            this.tabPage2.Controls.Add(this.label101);
            this.tabPage2.Controls.Add(this.label102);
            this.tabPage2.Controls.Add(this.dailyT_min);
            this.tabPage2.Controls.Add(this.dailyT_max);
            this.tabPage2.Controls.Add(this.dailyT_avg);
            this.tabPage2.Controls.Add(this.label97);
            this.tabPage2.Controls.Add(this.daily_n);
            this.tabPage2.Controls.Add(this.label96);
            this.tabPage2.Controls.Add(this.label93);
            this.tabPage2.Controls.Add(this.label89);
            this.tabPage2.Controls.Add(this.label40);
            this.tabPage2.Controls.Add(this.dailyET0);
            this.tabPage2.Controls.Add(this.dailyD);
            this.tabPage2.Controls.Add(this.dailyP);
            this.tabPage2.Location = new System.Drawing.Point(4, 22);
            this.tabPage2.Name = "tabPage2";
            this.tabPage2.Size = new System.Drawing.Size(803, 293);
            this.tabPage2.TabIndex = 10;
            this.tabPage2.Text = "Hydrological parameters";
            this.tabPage2.UseVisualStyleBackColor = true;
            // 
            // check_scaling_daily_weather
            // 
            this.check_scaling_daily_weather.AutoSize = true;
            this.check_scaling_daily_weather.Location = new System.Drawing.Point(125, 227);
            this.check_scaling_daily_weather.Name = "check_scaling_daily_weather";
            this.check_scaling_daily_weather.Size = new System.Drawing.Size(230, 17);
            this.check_scaling_daily_weather.TabIndex = 71;
            this.check_scaling_daily_weather.Text = "Scale daily weather with annual timeseries?";
            this.check_scaling_daily_weather.UseVisualStyleBackColor = true;
            // 
            // label106
            // 
            this.label106.AutoSize = true;
            this.label106.Location = new System.Drawing.Point(394, 114);
            this.label106.Name = "label106";
            this.label106.Size = new System.Drawing.Size(236, 13);
            this.label106.TabIndex = 70;
            this.label106.Text = "Snowfall and snowmelt temperature threshold [C]";
            // 
            // snow_threshold_textbox
            // 
            this.snow_threshold_textbox.Enabled = false;
            this.snow_threshold_textbox.Location = new System.Drawing.Point(340, 111);
            this.snow_threshold_textbox.Name = "snow_threshold_textbox";
            this.snow_threshold_textbox.Size = new System.Drawing.Size(40, 20);
            this.snow_threshold_textbox.TabIndex = 69;
            this.snow_threshold_textbox.Text = "0";
            this.snow_threshold_textbox.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // label105
            // 
            this.label105.AutoSize = true;
            this.label105.Location = new System.Drawing.Point(394, 76);
            this.label105.Name = "label105";
            this.label105.Size = new System.Drawing.Size(174, 13);
            this.label105.TabIndex = 68;
            this.label105.Text = "Snowmelt factor [m degree-1 day-1]";
            // 
            // snowmelt_factor_textbox
            // 
            this.snowmelt_factor_textbox.Enabled = false;
            this.snowmelt_factor_textbox.Location = new System.Drawing.Point(340, 73);
            this.snowmelt_factor_textbox.Name = "snowmelt_factor_textbox";
            this.snowmelt_factor_textbox.Size = new System.Drawing.Size(40, 20);
            this.snowmelt_factor_textbox.TabIndex = 67;
            this.snowmelt_factor_textbox.Text = "0.004";
            this.snowmelt_factor_textbox.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // label104
            // 
            this.label104.AutoSize = true;
            this.label104.Location = new System.Drawing.Point(333, 15);
            this.label104.Name = "label104";
            this.label104.Size = new System.Drawing.Size(118, 13);
            this.label104.TabIndex = 66;
            this.label104.Text = "Properties of study area";
            // 
            // latitude_min
            // 
            this.latitude_min.Enabled = false;
            this.latitude_min.Location = new System.Drawing.Point(397, 35);
            this.latitude_min.Name = "latitude_min";
            this.latitude_min.Size = new System.Drawing.Size(44, 20);
            this.latitude_min.TabIndex = 65;
            this.latitude_min.Text = "22";
            this.latitude_min.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // label103
            // 
            this.label103.AutoSize = true;
            this.label103.Location = new System.Drawing.Point(446, 38);
            this.label103.Name = "label103";
            this.label103.Size = new System.Drawing.Size(100, 13);
            this.label103.TabIndex = 64;
            this.label103.Text = "Latitude [deg], [min]";
            // 
            // latitude_deg
            // 
            this.latitude_deg.Enabled = false;
            this.latitude_deg.Location = new System.Drawing.Point(340, 35);
            this.latitude_deg.Name = "latitude_deg";
            this.latitude_deg.Size = new System.Drawing.Size(40, 20);
            this.latitude_deg.TabIndex = 63;
            this.latitude_deg.Text = "53";
            this.latitude_deg.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // label100
            // 
            this.label100.AutoSize = true;
            this.label100.Location = new System.Drawing.Point(143, 148);
            this.label100.Name = "label100";
            this.label100.Size = new System.Drawing.Size(59, 13);
            this.label100.TabIndex = 62;
            this.label100.Text = "Daily T min";
            // 
            // label101
            // 
            this.label101.AutoSize = true;
            this.label101.Location = new System.Drawing.Point(143, 174);
            this.label101.Name = "label101";
            this.label101.Size = new System.Drawing.Size(62, 13);
            this.label101.TabIndex = 61;
            this.label101.Text = "Daily T max";
            // 
            // label102
            // 
            this.label102.AutoSize = true;
            this.label102.Location = new System.Drawing.Point(143, 117);
            this.label102.Name = "label102";
            this.label102.Size = new System.Drawing.Size(82, 13);
            this.label102.TabIndex = 60;
            this.label102.Text = "Daily T average";
            // 
            // dailyT_min
            // 
            this.dailyT_min.Enabled = false;
            this.dailyT_min.Location = new System.Drawing.Point(37, 145);
            this.dailyT_min.Name = "dailyT_min";
            this.dailyT_min.Size = new System.Drawing.Size(100, 20);
            this.dailyT_min.TabIndex = 59;
            this.dailyT_min.Text = "D:\\PhD\\projects\\1g_basic LORICA development\\daily water\\Grunow\\Tminday_grunow.csv" +
    "";
            this.dailyT_min.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // dailyT_max
            // 
            this.dailyT_max.Enabled = false;
            this.dailyT_max.Location = new System.Drawing.Point(37, 171);
            this.dailyT_max.Name = "dailyT_max";
            this.dailyT_max.Size = new System.Drawing.Size(100, 20);
            this.dailyT_max.TabIndex = 58;
            this.dailyT_max.Text = "D:\\PhD\\projects\\1g_basic LORICA development\\daily water\\Grunow\\Tmaxday_grunow.csv" +
    "";
            this.dailyT_max.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // dailyT_avg
            // 
            this.dailyT_avg.Enabled = false;
            this.dailyT_avg.Location = new System.Drawing.Point(37, 114);
            this.dailyT_avg.Name = "dailyT_avg";
            this.dailyT_avg.Size = new System.Drawing.Size(100, 20);
            this.dailyT_avg.TabIndex = 57;
            this.dailyT_avg.Text = "D:\\PhD\\projects\\1g_basic LORICA development\\daily water\\Grunow\\Tavgday_grunow.csv" +
    "";
            this.dailyT_avg.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // label97
            // 
            this.label97.AutoSize = true;
            this.label97.Location = new System.Drawing.Point(142, 204);
            this.label97.Name = "label97";
            this.label97.Size = new System.Drawing.Size(83, 13);
            this.label97.TabIndex = 56;
            this.label97.Text = "Amount of years";
            // 
            // daily_n
            // 
            this.daily_n.Enabled = false;
            this.daily_n.Location = new System.Drawing.Point(36, 201);
            this.daily_n.Name = "daily_n";
            this.daily_n.Size = new System.Drawing.Size(100, 20);
            this.daily_n.TabIndex = 55;
            this.daily_n.Text = "6";
            this.daily_n.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // label96
            // 
            this.label96.AutoSize = true;
            this.label96.Location = new System.Drawing.Point(20, 15);
            this.label96.Name = "label96";
            this.label96.Size = new System.Drawing.Size(229, 13);
            this.label96.TabIndex = 54;
            this.label96.Text = "Insert a series of yearly records of the following:";
            // 
            // label93
            // 
            this.label93.AutoSize = true;
            this.label93.Location = new System.Drawing.Point(142, 65);
            this.label93.Name = "label93";
            this.label93.Size = new System.Drawing.Size(53, 13);
            this.label93.TabIndex = 53;
            this.label93.Text = "Daily ET0";
            // 
            // label89
            // 
            this.label89.AutoSize = true;
            this.label89.Location = new System.Drawing.Point(142, 91);
            this.label89.Name = "label89";
            this.label89.Size = new System.Drawing.Size(71, 13);
            this.label89.TabIndex = 52;
            this.label89.Text = "Daily duration";
            // 
            // label40
            // 
            this.label40.AutoSize = true;
            this.label40.Location = new System.Drawing.Point(142, 34);
            this.label40.Name = "label40";
            this.label40.Size = new System.Drawing.Size(40, 13);
            this.label40.TabIndex = 51;
            this.label40.Text = "Daily P";
            // 
            // dailyET0
            // 
            this.dailyET0.Enabled = false;
            this.dailyET0.Location = new System.Drawing.Point(36, 62);
            this.dailyET0.Name = "dailyET0";
            this.dailyET0.Size = new System.Drawing.Size(100, 20);
            this.dailyET0.TabIndex = 50;
            this.dailyET0.Text = "D:\\PhD\\projects\\1g_basic LORICA development\\daily water\\Grunow\\ET0day_grunow.csv";
            this.dailyET0.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // dailyD
            // 
            this.dailyD.Enabled = false;
            this.dailyD.Location = new System.Drawing.Point(36, 88);
            this.dailyD.Name = "dailyD";
            this.dailyD.Size = new System.Drawing.Size(100, 20);
            this.dailyD.TabIndex = 49;
            this.dailyD.Text = "D:\\PhD\\projects\\1g_basic LORICA development\\daily water\\Grunow\\Dday_grunow.csv";
            this.dailyD.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // dailyP
            // 
            this.dailyP.Enabled = false;
            this.dailyP.Location = new System.Drawing.Point(36, 31);
            this.dailyP.Name = "dailyP";
            this.dailyP.Size = new System.Drawing.Size(100, 20);
            this.dailyP.TabIndex = 48;
            this.dailyP.Text = "D:\\PhD\\projects\\1g_basic LORICA development\\daily water\\Grunow\\Pday_grunow.csv";
            this.dailyP.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // view_maps_checkbox
            // 
            this.view_maps_checkbox.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Left)));
            this.view_maps_checkbox.AutoSize = true;
            this.view_maps_checkbox.Checked = true;
            this.view_maps_checkbox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.view_maps_checkbox.Location = new System.Drawing.Point(354, 398);
            this.view_maps_checkbox.Name = "view_maps_checkbox";
            this.view_maps_checkbox.Size = new System.Drawing.Size(86, 17);
            this.view_maps_checkbox.TabIndex = 155;
            this.view_maps_checkbox.Text = "make maps?";
            this.view_maps_checkbox.UseVisualStyleBackColor = true;
            // 
            // Mother_form
            // 
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.None;
            this.AutoScroll = true;
            this.AutoSizeMode = System.Windows.Forms.AutoSizeMode.GrowAndShrink;
            this.ClientSize = new System.Drawing.Size(1184, 497);
            this.Controls.Add(this.view_maps_checkbox);
            this.Controls.Add(this.map_controls);
            this.Controls.Add(this.View_tabs_checkbox);
            this.Controls.Add(this.End_button);
            this.Controls.Add(this.start_button);
            this.Controls.Add(this.tabControl1);
            this.Controls.Add(this.statusBar1);
            this.Icon = ((System.Drawing.Icon)(resources.GetObject("$this.Icon")));
            this.MaximumSize = new System.Drawing.Size(1200, 700);
            this.Menu = this.mainMenu1;
            this.MinimumSize = new System.Drawing.Size(823, 530);
            this.Name = "Mother_form";
            this.StartPosition = System.Windows.Forms.FormStartPosition.CenterScreen;
            this.Text = "Timer";
            this.Load += new System.EventHandler(this.Form1_Load);
            Landsliding.ResumeLayout(false);
            Landsliding.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox4)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.InfoStatusPanel)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.TimeStatusPanel)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.ProcessStatusPanel)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.out_sed_statuspanel)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.total_tillage_statuspanel)).EndInit();
            this.groupBox13.ResumeLayout(false);
            this.groupBox13.PerformLayout();
            this.groupBox3.ResumeLayout(false);
            this.groupBox9.ResumeLayout(false);
            this.groupBox9.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.trackBar1)).EndInit();
            this.map_controls.ResumeLayout(false);
            this.map_controls.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.trackBar2)).EndInit();
            this.Output.ResumeLayout(false);
            this.groupBox6.ResumeLayout(false);
            this.groupBox12.ResumeLayout(false);
            this.groupBox12.PerformLayout();
            this.groupBox11.ResumeLayout(false);
            this.groupBox11.PerformLayout();
            this.groupBox1.ResumeLayout(false);
            this.groupBox1.PerformLayout();
            this.groupBox5.ResumeLayout(false);
            this.groupBox5.PerformLayout();
            this.UTMgroupBox.ResumeLayout(false);
            this.UTMgroupBox.PerformLayout();
            this.Run.ResumeLayout(false);
            this.Run.PerformLayout();
            this.groupBox2.ResumeLayout(false);
            this.groupBox2.PerformLayout();
            this.groupBox7.ResumeLayout(false);
            this.groupBox7.PerformLayout();
            this.Input.ResumeLayout(false);
            this.Input.PerformLayout();
            this.groupBox8.ResumeLayout(false);
            this.groupBox8.PerformLayout();
            this.Processes.ResumeLayout(false);
            this.Process_tabs.ResumeLayout(false);
            this.Water.ResumeLayout(false);
            this.Water.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox1)).EndInit();
            this.Tillage.ResumeLayout(false);
            this.Tillage.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox2)).EndInit();
            this.Creeper.ResumeLayout(false);
            this.Creeper.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox3)).EndInit();
            this.Solifluction.ResumeLayout(false);
            this.Solifluction.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox5)).EndInit();
            this.Rock_weathering.ResumeLayout(false);
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox6)).EndInit();
            this.groupBox10.ResumeLayout(false);
            this.groupBox10.PerformLayout();
            this.Tectonics.ResumeLayout(false);
            this.groupBox14.ResumeLayout(false);
            this.groupBox14.PerformLayout();
            this.groupBox16.ResumeLayout(false);
            this.groupBox16.PerformLayout();
            this.groupBox4.ResumeLayout(false);
            this.groupBox4.PerformLayout();
            this.groupBox15.ResumeLayout(false);
            this.groupBox15.PerformLayout();
            this.treefall.ResumeLayout(false);
            this.treefall.PerformLayout();
            this.tabControl1.ResumeLayout(false);
            this.tabPage1.ResumeLayout(false);
            this.tabControl2.ResumeLayout(false);
            this.physical.ResumeLayout(false);
            this.physical.PerformLayout();
            this.chemical.ResumeLayout(false);
            this.chemical.PerformLayout();
            this.clay.ResumeLayout(false);
            this.clay.PerformLayout();
            this.bioturbation.ResumeLayout(false);
            this.bioturbation.PerformLayout();
            this.carbon.ResumeLayout(false);
            this.carbon.PerformLayout();
            this.decalcification.ResumeLayout(false);
            this.decalcification.PerformLayout();
            this.tabPage2.ResumeLayout(false);
            this.tabPage2.PerformLayout();
            this.ResumeLayout(false);
            this.PerformLayout();

        }
        #endregion

        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        static void Main()
        {
            Application.Run(new Mother_form());
        }   // creates the forms

        LORICA4.Output_timeseries timeseries = new LORICA4.Output_timeseries();
        LORICA4.Output_profile profile = new LORICA4.Output_profile();
        LORICA4.Landuse_determinator landuse_determinator = new LORICA4.Landuse_determinator();
        LORICA4.About_LORICA aboutbox = new LORICA4.About_LORICA();
        public LORICA4.Soil_specifier soildata = new LORICA4.Soil_specifier();

        #region memory, reading and writing utilities

        private void menuItemConfigFileOpen_Click(object sender, System.EventArgs e)
        {
            //opens a runfile
            XmlTextReader xreader = null;
            int read_error = 0;
            OpenFileDialog openFileDialog1 = new OpenFileDialog();

            openFileDialog1.InitialDirectory = GlobalMethods.Workdir;
            openFileDialog1.Filter = "cfg files (*.xml)|*.xml|All files (*.*)|*.*";
            openFileDialog1.FilterIndex = 1;
            openFileDialog1.RestoreDirectory = false;

            
            var t = new Thread((ThreadStart)(() => {
                if (openFileDialog1.ShowDialog() == DialogResult.OK)
                {
                    cfgname = openFileDialog1.FileName;

                    xreader = new XmlTextReader(cfgname);

                }
                else return;
                
            }));

            t.SetApartmentState(ApartmentState.STA);
            t.Start();
            t.Join();
            //Read the file
            if (xreader != null)
            {
                try { xreader.ReadStartElement("Parms"); }
                catch { read_error = 1; }
                try { xreader.ReadStartElement("Processes"); }
                catch { read_error = 1; }
                try { xreader.ReadStartElement("Water_erosion"); }
                catch { read_error = 1; }
                try { Water_ero_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("water_active")); }
                catch { read_error = 1; }
                try { parameter_m_textbox.Text = xreader.ReadElementString("para_m"); }
                catch { read_error = 1; }
                try { parameter_n_textbox.Text = xreader.ReadElementString("para_n"); }
                catch { read_error = 1; }
                try { parameter_conv_textbox.Text = xreader.ReadElementString("para_p"); }
                catch { read_error = 1; }
                try { parameter_K_textbox.Text = xreader.ReadElementString("para_K"); }
                catch { read_error = 1; }
                try { erosion_threshold_textbox.Text = xreader.ReadElementString("para_ero_threshold"); }
                catch { read_error = 1; }
                try { rock_protection_constant_textbox.Text = xreader.ReadElementString("para_rock_protection_const"); }
                catch { read_error = 1; }
                try { bio_protection_constant_textbox.Text = xreader.ReadElementString("para_bio_protection_const"); }
                catch { read_error = 1; }
                try { selectivity_constant_textbox.Text = xreader.ReadElementString("para_selectivity"); }
                catch { read_error = 1; }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Tillage"); }
                catch { read_error = 1; }
                try { Tillage_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("tillage_active")); }
                catch { read_error = 1; }
                try { parameter_ploughing_depth_textbox.Text = xreader.ReadElementString("para_plough_depth"); }
                catch { read_error = 1; }
                try { parameter_tillage_constant_textbox.Text = xreader.ReadElementString("para_tillage_constant"); }
                catch { read_error = 1; }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Weathering"); }
                catch { read_error = 1; }
                try { Biological_weathering_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("bio_weathering_active")); }
                catch { read_error = 1; }
                try { parameter_P0_textbox.Text = xreader.ReadElementString("para_P0"); }
                catch { read_error = 1; }
                try { parameter_k1_textbox.Text = xreader.ReadElementString("para_k1"); }
                catch { read_error = 1; }
                try { parameter_k2_textbox.Text = xreader.ReadElementString("para_k2"); }
                catch { read_error = 1; }
                try { parameter_Pa_textbox.Text = xreader.ReadElementString("para_Pa"); }
                catch { read_error = 1; }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Landsliding"); }
                catch { read_error = 1; }
                try { Landslide_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("landsliding_active")); }
                catch { read_error = 1; }
                try { radio_ls_absolute.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("radio_ls_absolute")); }
                catch { read_error = 1; }
                try { radio_ls_fraction.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("radio_ls_fraction")); }
                catch { read_error = 1; }
                try { text_ls_abs_rain_intens.Text = xreader.ReadElementString("para_absolute_rain_intens"); }
                catch { read_error = 1; }
                try { text_ls_rel_rain_intens.Text = xreader.ReadElementString("para_relative_rain_intens"); }
                catch { read_error = 1; }
                try { textBox_ls_coh.Text = xreader.ReadElementString("para_cohesion"); }
                catch { read_error = 1; }
                try { textBox_ls_ifr.Text = xreader.ReadElementString("para_friction"); }
                catch { read_error = 1; }
                try { textBox_ls_bd.Text = xreader.ReadElementString("para_density"); }
                catch { read_error = 1; }
                try { textBox_ls_trans.Text = xreader.ReadElementString("para_transmissivity"); }
                catch { read_error = 1; }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Creep"); }
                catch { read_error = 1; }
                try { creep_active_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("creep_active")); }
                catch { read_error = 1; }
                try { parameter_diffusivity_textbox.Text = xreader.ReadElementString("para_diffusivity"); }
                catch { read_error = 1; }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Tree_fall"); }
                catch { read_error = 1; }
                try { treefall_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("treefall_active")); }
                catch { read_error = 1; }
                try { tf_W.Text = xreader.ReadElementString("tf_width"); }
                catch { read_error = 1; }
                try { tf_D.Text = xreader.ReadElementString("tf_depth"); }
                catch { read_error = 1; }
                try { tf_growth.Text = xreader.ReadElementString("tf_growth"); }
                catch { read_error = 1; }
                try { tf_age.Text = xreader.ReadElementString("tf_age"); }
                catch { read_error = 1; }
                try { tf_freq.Text = xreader.ReadElementString("tf_freq"); }
                catch { read_error = 1; }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }


                try { xreader.ReadStartElement("Soil_forming_processes"); }
                catch { read_error = 1; }
                try { xreader.ReadStartElement("Physical_weathering"); }
                catch { read_error = 1; }
                try { soil_phys_weath_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("phys_weath_active")); }
                catch { read_error = 1; }
                try { Physical_weath_C1_textbox.Text = xreader.ReadElementString("weath_rate_constant"); }
                catch { read_error = 1; }
                try { physical_weath_constant1.Text = xreader.ReadElementString("constant1"); }
                catch { read_error = 1; }
                try { physical_weath_constant2.Text = xreader.ReadElementString("constant2"); }
                catch { read_error = 1; }
                try { upper_particle_coarse_textbox.Text = xreader.ReadElementString("size_coarse"); }
                catch { read_error = 1; }
                try { upper_particle_sand_textbox.Text = xreader.ReadElementString("size_sand"); }
                catch { read_error = 1; }
                try { upper_particle_silt_textbox.Text = xreader.ReadElementString("size_silt"); }
                catch { read_error = 1; }
                try { upper_particle_clay_textbox.Text = xreader.ReadElementString("size_clay"); }
                catch { read_error = 1; }
                try { upper_particle_fine_clay_textbox.Text = xreader.ReadElementString("size_fine"); }
                catch { read_error = 1; }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Chemical_weathering"); }
                catch { read_error = 1; }
                try { soil_chem_weath_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("chem_weath_active")); }
                catch { read_error = 1; }
                try { chem_weath_rate_constant_textbox.Text = xreader.ReadElementString("weath_rate_constant"); }
                catch { read_error = 1; }
                try { chem_weath_depth_constant_textbox.Text = xreader.ReadElementString("constant3"); }
                catch { read_error = 1; }
                try { chem_weath_specific_coefficient_textbox.Text = xreader.ReadElementString("constant4"); }
                catch { read_error = 1; }
                try { specific_area_coarse_textbox.Text = xreader.ReadElementString("surface_coarse"); }
                catch { read_error = 1; }
                try { specific_area_sand_textbox.Text = xreader.ReadElementString("surface_sand"); }
                catch { read_error = 1; }
                try { specific_area_silt_textbox.Text = xreader.ReadElementString("surface_silt"); }
                catch { read_error = 1; }
                try { specific_area_clay_textbox.Text = xreader.ReadElementString("surface_clay"); }
                catch { read_error = 1; }
                try { specific_area_fine_clay_textbox.Text = xreader.ReadElementString("surface_fine_clay"); }
                catch { read_error = 1; }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Clay_dynamics"); }
                catch { read_error = 1; }
                try { soil_clay_transloc_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("clay_dynamics_active")); }
                catch { read_error = 1; }
                try { clay_neoform_constant_textbox.Text = xreader.ReadElementString("neoform_rate_constant"); }
                catch { read_error = 1; }
                try { clay_neoform_C1_textbox.Text = xreader.ReadElementString("constant5"); }
                catch { read_error = 1; }
                try { clay_neoform_C2_textbox.Text = xreader.ReadElementString("constant6"); }
                catch { read_error = 1; }
                try { maximum_eluviation_textbox.Text = xreader.ReadElementString("max_eluviation"); }
                catch { read_error = 1; }
                try { eluviation_coefficient_textbox.Text = xreader.ReadElementString("eluviation_coefficient"); }
                catch { read_error = 1; Debug.WriteLine("xml1"); }
                try { ct_Jagercikova.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("ct_Jagercikova_active")); } //MMxml
                catch { read_error = 1; }
                try { ct_v0_Jagercikova.Text = xreader.ReadElementString("ct_v0_Jagercikova"); } //MMxml
                catch { read_error = 2; Debug.WriteLine("xml2"); }
                try { ct_dd_Jagercikova.Text = xreader.ReadElementString("ct_dd_Jagercikova"); } //MMxml
                catch { read_error = 2; Debug.WriteLine("xml3"); }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Bioturbation"); }
                catch { read_error = 1; }
                try { soil_bioturb_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("bioturbation_active")); }
                catch { read_error = 1; }
                try { potential_bioturbation_textbox.Text = xreader.ReadElementString("potential_bioturb"); }
                catch { read_error = 1; }
                try { bioturbation_depth_decay_textbox.Text = xreader.ReadElementString("bioturb_depth_decay"); }
                catch { read_error = 1; }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Carboncycle"); }
                catch { read_error = 1; }
                try { soil_carbon_cycle_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("carboncycle_active")); }
                catch { read_error = 1; }
                try { carbon_input_textbox.Text = xreader.ReadElementString("carbon_input"); }
                catch { read_error = 1; }
                try { carbon_depth_decay_textbox.Text = xreader.ReadElementString("carbon_depth_decay"); }
                catch { read_error = 1; }
                try { carbon_humification_fraction_textbox.Text = xreader.ReadElementString("carbon_hum_fraction"); }
                catch { read_error = 1; }
                try { carbon_y_decomp_rate_textbox.Text = xreader.ReadElementString("carbon_y_decomp"); }
                catch { read_error = 1; }
                try { carbon_y_depth_decay_textbox.Text = xreader.ReadElementString("carbon_y_depth_decay"); }
                catch { read_error = 1; }
                try { carbon_y_twi_decay_textbox.Text = xreader.ReadElementString("carbon_y_twi_decay"); }
                catch { read_error = 1; }
                try { carbon_o_decomp_rate_textbox.Text = xreader.ReadElementString("carbon_o_decomp"); }
                catch { read_error = 1; }
                try { carbon_o_depth_decay_textbox.Text = xreader.ReadElementString("carbon_o_depth_decay"); }
                catch { read_error = 1; }
                try { carbon_o_twi_decay_textbox.Text = xreader.ReadElementString("carbon_o_twi_decay"); }
                catch { read_error = 1; }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Inputs"); }
                catch { read_error = 1; }
                try { check_space_DTM.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("check_space_DTM")); }
                catch { read_error = 1; }
                try { check_space_soildepth.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("check_space_soil")); }
                catch { read_error = 1; }
                try { check_space_landuse.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("check_space_landuse")); }
                catch { read_error = 1; }
                try { check_space_till_fields.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("check_space_tillfields")); }
                catch { read_error = 1; }
                try { check_space_rain.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("check_space_rain")); }
                catch { read_error = 1; }
                try { check_space_infil.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("check_space_infil")); }
                catch { read_error = 1; }
                try { check_space_evap.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("check_space_evap")); }
                catch { read_error = 1; }
                try { check_time_landuse.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("check_time_landuse")); }
                catch { read_error = 1; }
                try { check_time_till_fields.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("check_time_tillfields")); }
                catch { read_error = 1; }
                try { check_time_rain.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("check_time_rain")); }
                catch { read_error = 1; }
                try { check_time_infil.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("check_time_infil")); }
                catch { read_error = 1; }
                try { check_time_evap.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("check_time_evap")); }
                catch { read_error = 1; }

                try { daily_water.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("dailywater")); } //MMxml
                catch { read_error = 1; }
                try { dailyP.Text = xreader.ReadElementString("dailyP"); }//MMxml
                catch { read_error = 1; Debug.WriteLine("xml4"); }
                try { dailyET0.Text = xreader.ReadElementString("dailyET0"); }//MMxml
                catch { read_error = 1; }
                try { dailyD.Text = xreader.ReadElementString("dailyD"); }//MMxml
                catch { read_error = 1; }
                try { dailyT_avg.Text = xreader.ReadElementString("dailyT_avg"); }//MMxml
                catch { read_error = 1; }
                try { dailyT_min.Text = xreader.ReadElementString("dailyT_min"); }//MMxml
                catch { read_error = 1; Debug.WriteLine("xml5"); }
                try { dailyT_max.Text = xreader.ReadElementString("dailyT_max"); }//MMxml
                catch { read_error = 1; }
                try { latitude_deg.Text = xreader.ReadElementString("latitude_deg"); }//MMxml
                catch { read_error = 1; }
                try { latitude_min.Text = xreader.ReadElementString("latitude_min"); }//MMxml
                catch { read_error = 1; Debug.WriteLine("xml6"); }
                try { snowmelt_factor_textbox.Text = xreader.ReadElementString("snowmelt_factor"); }//MMxml
                catch { read_error = 1; }
                try { snow_threshold_textbox.Text = xreader.ReadElementString("snowmelt_threshold"); }//MMxml
                catch { read_error = 1; Debug.WriteLine("xml7"); }
                try { daily_n.Text = xreader.ReadElementString("daily_n_years"); }//MMxml
                catch { read_error = 1; Debug.WriteLine("xml7.1"); }
                try { check_scaling_daily_weather.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("scaledailyweather")); } //MMxml
                catch { read_error = 1; Debug.WriteLine("xml7.2"); }



                try { dtm_input_filename_textbox.Text = xreader.ReadElementString("dtm_input_filename"); }
                catch { read_error = 1; }
                try { soildepth_input_filename_textbox.Text = xreader.ReadElementString("soildepth_input_filename"); }
                catch { read_error = 1; }
                try { landuse_input_filename_textbox.Text = xreader.ReadElementString("landuse_input_filename"); }
                catch { read_error = 1; }
                try { tillfields_input_filename_textbox.Text = xreader.ReadElementString("tillfields_input_filename"); }
                catch { read_error = 1; }
                try { rain_input_filename_textbox.Text = xreader.ReadElementString("rain_input_filename"); }
                catch { read_error = 1; }
                try { infil_input_filename_textbox.Text = xreader.ReadElementString("infil_input_filename"); }
                catch { read_error = 1; Debug.WriteLine("xml8"); }
                try { evap_input_filename_textbox.Text = xreader.ReadElementString("evap_input_filename"); }
                catch { read_error = 1; }

                try { soildepth_constant_value_box.Text = xreader.ReadElementString("soildepth_constant_value"); }
                catch { read_error = 1; }
                try { landuse_constant_value_box.Text = xreader.ReadElementString("landuse_constant_value"); }
                catch { read_error = 1; }
                try { tillfields_constant_textbox.Text = xreader.ReadElementString("tillfields_constant_value"); }
                catch { read_error = 1; }
                try { rainfall_constant_value_box.Text = xreader.ReadElementString("rain_constant_value"); }
                catch { read_error = 1; }
                try { infil_constant_value_box.Text = xreader.ReadElementString("infil_constant_value"); }
                catch { read_error = 1; }
                try { evap_constant_value_box.Text = xreader.ReadElementString("evap_constant_value"); }
                catch { read_error = 1; }

                try { fill_sinks_before_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("check_fill_sinks_before")); }
                catch { read_error = 1; Debug.WriteLine("xml9"); }
                try { fill_sinks_during_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("check_fill_sinks_during")); }
                catch { read_error = 1; }

                try { xreader.ReadEndElement(); }
                catch { read_error = 1; Debug.WriteLine("xm20"); }

                try { xreader.ReadStartElement("Run"); }
                catch { read_error = 1; Debug.WriteLine("xm21"); }
                try { runs_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("runs_radiobutton")); }
                catch { read_error = 1; Debug.WriteLine("xm22"); }
                try { Number_runs_textbox.Text = xreader.ReadElementString("number_runs"); }
                catch { read_error = 1; Debug.WriteLine("xm23"); }


                try { xreader.ReadStartElement("Specialsettings"); }
                catch { read_error = 1; Debug.WriteLine("xm24"); }
                try { Ik_ben_Marijn.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("Spitsbergen")); }
                catch { read_error = 1; Debug.WriteLine("xm25"); }
                try { version_lux_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("Luxembourg")); }
                catch { read_error = 1; Debug.WriteLine("xm26"); }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; Debug.WriteLine("xm27"); }

                try { xreader.ReadStartElement("CalibrationSensitivity"); }
                catch { read_error = 2; }
                try
                {
                    Calibration_button.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("calibration_active_button"));
                    calibration_ratios_textbox.Text = xreader.ReadElementString("calibration_ratios_string");
                    calibration_levels_textbox.Text = xreader.ReadElementString("calibration_levels");
                    calibration_ratio_reduction_parameter_textbox.Text = xreader.ReadElementString("calibration_ratio_reduction_per_level");
                    xreader.ReadEndElement();
                }
                catch { read_error = 2; Debug.WriteLine("xm28"); }

                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Output"); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("File_Output"); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Moment_of_Output"); }
                catch { read_error = 1; }
                try { Final_output_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("final_output_checkbox")); }
                catch { read_error = 1; }
                try { Regular_output_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("regular_output_checkbox")); }
                catch { read_error = 1; }
                try { Box_years_output.Text = xreader.ReadElementString("years_between_outputs"); }
                catch { read_error = 1; }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Type_of_Output"); }
                catch { read_error = 1; }
                try { cumulative_output_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("cumulative")); }
                catch { read_error = 1; }
                try { annual_output_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("annual")); }
                catch { read_error = 1; }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Maps_to_Output"); }
                catch { read_error = 1; }
                try { Altitude_output_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("alti")); }
                catch { read_error = 1; }
                try { Alt_change_output_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("altichange")); }
                catch { read_error = 1; }
                try { Soildepth_output_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("soildepth")); }
                catch { read_error = 1; }
                try { all_process_output_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("all_processes")); }
                catch { read_error = 1; }
                try { water_output_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("waterflow")); }
                catch { read_error = 1; }
                try { depressions_output_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("depressions")); }
                catch { read_error = 1; }
                try { diagnostic_output_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("diagnostics")); }
                catch { read_error = 1; }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Other_outputs"); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Avi_Animation"); }
                catch { read_error = 1; }
                try { checkBoxGenerateAVIFile.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("Save_as_AVI")); }
                catch { read_error = 1; }
                try { textBoxAVIFile.Text = xreader.ReadElementString("AVIfilename"); }
                catch { read_error = 1; }
                try { saveintervalbox.Text = xreader.ReadElementString("years_between_AVI_outputs"); }
                catch { read_error = 1; }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Google_Animation"); }
                catch { read_error = 1; }
                try { googleAnimationCheckbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("Save_as_Google")); }
                catch { read_error = 1; }
                try { googleAnimationTextBox.Text = xreader.ReadElementString("Google_filename"); }
                catch { read_error = 1; }
                try { googAnimationSaveInterval.Text = xreader.ReadElementString("years_between_Google_outputs"); }
                catch { read_error = 1; }
                try { googleBeginDate.Text = xreader.ReadElementString("begindate"); }
                catch { read_error = 1; }
                try { UTMgridcheckbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("UTM")); }
                catch { read_error = 1; }
                try { UTMsouthcheck.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("South")); }
                catch { read_error = 1; }
                try { UTMzonebox.Text = xreader.ReadElementString("UTMzone"); }
                catch { read_error = 1; }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Timeseries"); }
                catch { read_error = 1; }
                try { timeseries.timeseries_total_ero_check.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("total_erosion")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_total_dep_check.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("total_deposition")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_net_ero_check.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("net_erosion")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_SDR_check.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("SDR")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_total_average_alt_check.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("total_average_alt")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_total_rain_check.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("total_rain")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_total_infil_check.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("total_infil")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_total_evap_check.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("total_evap")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_total_outflow_check.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("total_outflow")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_number_waterflow_check.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("wet_cells")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_number_erosion_check.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("eroded_cells")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_number_dep_check.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("deposited_cells")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_cell_altitude_check.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("cell_altitude")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_cell_waterflow_check.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("cell_waterflow")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_textbox_waterflow_threshold.Text = xreader.ReadElementString("waterflow_threshold"); }
                catch { read_error = 1; }
                try { timeseries.timeseries_textbox_erosion_threshold.Text = xreader.ReadElementString("erosion_threshold"); }
                catch { read_error = 1; }
                try { timeseries.timeseries_textbox_deposition_threshold.Text = xreader.ReadElementString("deposition_threshold"); }
                catch { read_error = 1; }
                try { timeseries.timeseries_textbox_cell_row.Text = xreader.ReadElementString("cell_row"); }
                catch { read_error = 1; }
                try { timeseries.timeseries_textbox_cell_col.Text = xreader.ReadElementString("cell_col"); }
                catch { read_error = 1; }



                try { timeseries.total_OM_input_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("total_OM_input")); }
                catch { read_error = 1; }
                try { timeseries.total_average_soilthickness_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("total_average_soil_thickness")); }
                catch { read_error = 1; }
                try { timeseries.total_phys_weath_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("total_phys_weath")); }
                catch { read_error = 1; }
                try { timeseries.total_chem_weath_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("total_chem_weath")); }
                catch { read_error = 1; }
                try { timeseries.total_fine_formed_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("total_fine_formed")); }
                catch { read_error = 1; }
                try { timeseries.total_fine_eluviated_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("total_fine_eluviated")); }
                catch { read_error = 1; }
                try { timeseries.total_mass_bioturbed_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("total_mass_bioturbed")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_soil_depth_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("timeseries_soil_depth")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_soil_mass_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("timeseries_soil_mass")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_coarser_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("timeseries_coarser")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_number_soil_thicker_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("timeseries_thicker")); }
                catch { read_error = 1; }
                try { timeseries.timeseries_soil_cell_col.Text = xreader.ReadElementString("soil_cell"); }
                catch { read_error = 1; }
                try { timeseries.timeseries_soil_cell_row.Text = xreader.ReadElementString("soil_col"); }
                catch { read_error = 1; }
                try { timeseries.timeseries_soil_coarser_fraction_textbox.Text = xreader.ReadElementString("coarser_fraction"); }
                catch { read_error = 1; }
                try { timeseries.timeseries_soil_thicker_textbox.Text = xreader.ReadElementString("thicker_threshold"); }
                catch { read_error = 1; }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Profiles"); }
                catch { read_error = 1; }
                try { profile.radio_pro1_row.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("profile1_row")); }
                catch { read_error = 1; }
                try { profile.radio_pro1_col.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("profile1_col")); }
                catch { read_error = 1; }
                try { profile.radio_pro2_row.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("profile2_row")); }
                catch { read_error = 1; }
                try { profile.radio_pro2_col.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("profile2_col")); }
                catch { read_error = 1; }
                try { profile.radio_pro3_row.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("profile3_row")); }
                catch { read_error = 1; }
                try { profile.radio_pro3_col.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("profile3_col")); }
                catch { read_error = 1; }
                try { profile.p1_row_col_box.Text = xreader.ReadElementString("p1_number"); }
                catch { read_error = 1; }
                try { profile.p2_row_col_box.Text = xreader.ReadElementString("p2_number"); }
                catch { read_error = 1; }
                try { profile.p3_row_col_box.Text = xreader.ReadElementString("p3_number"); }
                catch { read_error = 1; }
                try { profile.check_waterflow_profile1.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("p1_waterflow")); }
                catch { read_error = 1; }
                try { profile.check_altitude_profile1.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("p1_altitude")); }
                catch { read_error = 1; }
                try { profile.check_waterflow_profile2.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("p2_waterflow")); }
                catch { read_error = 1; }
                try { profile.check_altitude_profile2.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("p2_altitude")); }
                catch { read_error = 1; }
                try { profile.check_waterflow_profile3.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("p3_waterflow")); }
                catch { read_error = 1; }
                try { profile.check_altitude_profile3.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("p3_altitude")); }
                catch { read_error = 1; }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { view_maps_checkbox.Checked = XmlConvert.ToBoolean(xreader.ReadElementString("maps_required")); }
                catch { read_error = 1; }

                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                try { xreader.ReadStartElement("Soilfractions"); }
                catch { read_error = 1; }
                try { soildata.coarsebox.Text = xreader.ReadElementString("coarsefrac"); }
                catch { read_error = 1; }
                try { soildata.sandbox.Text = xreader.ReadElementString("sandfrac"); }
                catch { read_error = 1; }
                try { soildata.siltbox.Text = xreader.ReadElementString("siltfrac"); }
                catch { read_error = 1; }
                try { soildata.claybox.Text = xreader.ReadElementString("clayfrac"); }
                catch { read_error = 1; }
                try { soildata.fineclaybox.Text = xreader.ReadElementString("fclayfrac"); }
                catch { read_error = 1; }
                try { xreader.ReadEndElement(); }
                catch { read_error = 1; }

                if (read_error == 1) { MessageBox.Show("warning : not all runfile data could be read.\r\n LORICA can continue"); }
                if (read_error == 2) { MessageBox.Show("Error in new XML lines"); }


                xreader.Close();

                this.Text = basetext + " (" + Path.GetFileName(cfgname) + ")";
                start_button.Enabled = true;
                graphicToGoogleEarthButton.Visible = false;
                tabControl1.Visible = true;

            }
        }

        private void menuItemConfigFileSave_Click(object sender, System.EventArgs e)
        {
            XmlTextWriter xwriter;

            if ((sender == menuItemConfigFileSaveAs) || (cfgname == null))
            {

                SaveFileDialog saveFileDialog1 = new SaveFileDialog();

                saveFileDialog1.InitialDirectory = GlobalMethods.Workdir;
                saveFileDialog1.Filter = "cfg files (*.xml)|*.xml|All files (*.*)|*.*";
                saveFileDialog1.FilterIndex = 1;
                saveFileDialog1.RestoreDirectory = false;

                if (saveFileDialog1.ShowDialog() == DialogResult.OK)
                {
                    cfgname = saveFileDialog1.FileName;
                }
            }
            if (cfgname != null)
            {

                //Create a new XmlTextWriter.
                xwriter = new XmlTextWriter(cfgname, System.Text.Encoding.UTF8);
                //Write the beginning of the document including the 
                //document declaration. Standalone is true. 
                //Use indentation for readability.
                xwriter.Formatting = Formatting.Indented;
                xwriter.Indentation = 4;

                xwriter.WriteStartDocument(true);

                //Write the beginning of the "data" element. This is 
                //the opening tag to our data 
                xwriter.WriteStartElement("Parms");
                xwriter.WriteStartElement("Processes");
                xwriter.WriteStartElement("Water_erosion");
                xwriter.WriteElementString("water_active", XmlConvert.ToString(Water_ero_checkbox.Checked));
                xwriter.WriteElementString("para_m", parameter_m_textbox.Text);
                xwriter.WriteElementString("para_n", parameter_n_textbox.Text);
                xwriter.WriteElementString("para_p", parameter_conv_textbox.Text);
                xwriter.WriteElementString("para_K", parameter_K_textbox.Text);
                xwriter.WriteElementString("para_ero_threshold", erosion_threshold_textbox.Text);
                xwriter.WriteElementString("para_rock_protection_const", rock_protection_constant_textbox.Text);
                xwriter.WriteElementString("para_bio_protection_const", bio_protection_constant_textbox.Text);
                xwriter.WriteElementString("para_selectivity", selectivity_constant_textbox.Text);
                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Tillage");
                xwriter.WriteElementString("tillage_active", XmlConvert.ToString(Tillage_checkbox.Checked));
                xwriter.WriteElementString("para_plough_depth", parameter_ploughing_depth_textbox.Text);
                xwriter.WriteElementString("para_tillage_constant", parameter_tillage_constant_textbox.Text);
                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Weathering");
                xwriter.WriteElementString("bio_weathering_active", XmlConvert.ToString(Biological_weathering_checkbox.Checked));
                xwriter.WriteElementString("para_P0", parameter_P0_textbox.Text);
                xwriter.WriteElementString("para_k1", parameter_k1_textbox.Text);
                xwriter.WriteElementString("para_k2", parameter_k2_textbox.Text);
                xwriter.WriteElementString("para_Pa", parameter_Pa_textbox.Text);
                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Landsliding");
                xwriter.WriteElementString("landsliding_active", XmlConvert.ToString(Landslide_checkbox.Checked));
                xwriter.WriteElementString("radio_ls_absolute", XmlConvert.ToString(radio_ls_absolute.Checked));
                xwriter.WriteElementString("radio_ls_fraction", XmlConvert.ToString(radio_ls_fraction.Checked));
                xwriter.WriteElementString("para_absolute_rain_intens", text_ls_abs_rain_intens.Text);
                xwriter.WriteElementString("para_relative_rain_intens", text_ls_rel_rain_intens.Text);
                xwriter.WriteElementString("para_cohesion", textBox_ls_coh.Text);
                xwriter.WriteElementString("para_friction", textBox_ls_ifr.Text);
                xwriter.WriteElementString("para_density", textBox_ls_bd.Text);
                xwriter.WriteElementString("para_transmissivity", textBox_ls_trans.Text);
                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Creep");
                xwriter.WriteElementString("creep_active", XmlConvert.ToString(creep_active_checkbox.Checked));
                xwriter.WriteElementString("para_diffusivity", parameter_diffusivity_textbox.Text);
                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Tree_fall");
                xwriter.WriteElementString("treefall_active", XmlConvert.ToString(treefall_checkbox.Checked));
                xwriter.WriteElementString("tf_width", tf_W.Text);
                xwriter.WriteElementString("tf_depth", tf_D.Text);
                xwriter.WriteElementString("tf_growth", tf_growth.Text);
                xwriter.WriteElementString("tf_age", tf_age.Text);
                xwriter.WriteElementString("tf_freq", tf_freq.Text);
                xwriter.WriteEndElement();
                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Soil_forming_processes");
                xwriter.WriteStartElement("Physical_weathering");
                xwriter.WriteElementString("phys_weath_active", XmlConvert.ToString(soil_phys_weath_checkbox.Checked));
                xwriter.WriteElementString("weath_rate_constant", Physical_weath_C1_textbox.Text);
                xwriter.WriteElementString("constant1", physical_weath_constant1.Text);
                xwriter.WriteElementString("constant2", physical_weath_constant2.Text);
                xwriter.WriteElementString("size_coarse", upper_particle_coarse_textbox.Text);
                xwriter.WriteElementString("size_sand", upper_particle_sand_textbox.Text);
                xwriter.WriteElementString("size_silt", upper_particle_silt_textbox.Text);
                xwriter.WriteElementString("size_clay", upper_particle_clay_textbox.Text);
                xwriter.WriteElementString("size_fine", upper_particle_fine_clay_textbox.Text);
                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Chemical_weathering");
                xwriter.WriteElementString("chem_weath_active", XmlConvert.ToString(soil_chem_weath_checkbox.Checked));
                xwriter.WriteElementString("weath_rate_constant", chem_weath_rate_constant_textbox.Text);
                xwriter.WriteElementString("constant3", chem_weath_depth_constant_textbox.Text);
                xwriter.WriteElementString("constant4", chem_weath_specific_coefficient_textbox.Text);
                xwriter.WriteElementString("surface_coarse", specific_area_coarse_textbox.Text);
                xwriter.WriteElementString("surface_sand", specific_area_sand_textbox.Text);
                xwriter.WriteElementString("surface_silt", specific_area_silt_textbox.Text);
                xwriter.WriteElementString("surface_clay", specific_area_clay_textbox.Text);
                xwriter.WriteElementString("surface_fine_clay", specific_area_fine_clay_textbox.Text);
                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Clay_dynamics");
                xwriter.WriteElementString("clay_dynamics_active", XmlConvert.ToString(soil_clay_transloc_checkbox.Checked));
                xwriter.WriteElementString("neoform_rate_constant", clay_neoform_constant_textbox.Text);
                xwriter.WriteElementString("constant5", clay_neoform_C1_textbox.Text);
                xwriter.WriteElementString("constant6", clay_neoform_C2_textbox.Text);
                xwriter.WriteElementString("max_eluviation", maximum_eluviation_textbox.Text);
                xwriter.WriteElementString("eluviation_coefficient", eluviation_coefficient_textbox.Text);
                xwriter.WriteElementString("ct_Jagercikova_active", XmlConvert.ToString(ct_Jagercikova.Checked));
                xwriter.WriteElementString("ct_v0_Jagercikova", ct_v0_Jagercikova.Text);
                xwriter.WriteElementString("ct_dd_Jagercikova", ct_dd_Jagercikova.Text);
                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Bioturbation");
                xwriter.WriteElementString("bioturbation_active", XmlConvert.ToString(soil_bioturb_checkbox.Checked));
                xwriter.WriteElementString("potential_bioturb", potential_bioturbation_textbox.Text);
                xwriter.WriteElementString("bioturb_depth_decay", bioturbation_depth_decay_textbox.Text);
                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Carboncycle");
                xwriter.WriteElementString("carboncycle_active", XmlConvert.ToString(soil_carbon_cycle_checkbox.Checked));
                xwriter.WriteElementString("carbon_input", carbon_input_textbox.Text);
                xwriter.WriteElementString("carbon_depth_decay", carbon_depth_decay_textbox.Text);
                xwriter.WriteElementString("carbon_hum_fraction", carbon_humification_fraction_textbox.Text);
                xwriter.WriteElementString("carbon_y_decomp", carbon_y_decomp_rate_textbox.Text);
                xwriter.WriteElementString("carbon_y_depth_decay", carbon_y_depth_decay_textbox.Text);
                xwriter.WriteElementString("carbon_y_twi_decay", carbon_y_twi_decay_textbox.Text);
                xwriter.WriteElementString("carbon_o_decomp", carbon_o_decomp_rate_textbox.Text);
                xwriter.WriteElementString("carbon_o_depth_decay", carbon_o_depth_decay_textbox.Text);
                xwriter.WriteElementString("carbon_o_twi_decay", carbon_o_twi_decay_textbox.Text);
                xwriter.WriteEndElement();

                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Inputs");
                xwriter.WriteElementString("check_space_DTM", XmlConvert.ToString(check_space_DTM.Checked));
                xwriter.WriteElementString("check_space_soil", XmlConvert.ToString(check_space_soildepth.Checked));
                xwriter.WriteElementString("check_space_landuse", XmlConvert.ToString(check_space_landuse.Checked));
                xwriter.WriteElementString("check_space_tillfields", XmlConvert.ToString(check_space_till_fields.Checked));
                xwriter.WriteElementString("check_space_rain", XmlConvert.ToString(check_space_rain.Checked));
                xwriter.WriteElementString("check_space_infil", XmlConvert.ToString(check_space_infil.Checked));
                xwriter.WriteElementString("check_space_evap", XmlConvert.ToString(check_space_evap.Checked));
                xwriter.WriteElementString("check_time_landuse", XmlConvert.ToString(check_time_landuse.Checked));
                xwriter.WriteElementString("check_time_tillfields", XmlConvert.ToString(check_time_till_fields.Checked));
                xwriter.WriteElementString("check_time_rain", XmlConvert.ToString(check_time_rain.Checked));
                xwriter.WriteElementString("check_time_infil", XmlConvert.ToString(check_time_infil.Checked));
                xwriter.WriteElementString("check_time_evap", XmlConvert.ToString(check_time_evap.Checked));

                xwriter.WriteElementString("dailywater", XmlConvert.ToString(daily_water.Checked));
                xwriter.WriteElementString("dailyP", dailyP.Text);
                xwriter.WriteElementString("dailyET0", dailyET0.Text);
                xwriter.WriteElementString("dailyD", dailyD.Text);
                xwriter.WriteElementString("dailyT_avg", dailyT_avg.Text);
                xwriter.WriteElementString("dailyT_min", dailyT_min.Text);
                xwriter.WriteElementString("dailyT_max", dailyT_max.Text);
                xwriter.WriteElementString("latitude_deg", latitude_deg.Text);
                xwriter.WriteElementString("latitude_min", latitude_min.Text);
                xwriter.WriteElementString("snowmelt_factor", snowmelt_factor_textbox.Text);
                xwriter.WriteElementString("snowmelt_threshold", snow_threshold_textbox.Text);
                xwriter.WriteElementString("daily_n_years", daily_n.Text);
                xwriter.WriteElementString("scaledailyweather", XmlConvert.ToString(check_scaling_daily_weather.Checked));


                xwriter.WriteElementString("dtm_input_filename", dtm_input_filename_textbox.Text);
                xwriter.WriteElementString("soildepth_input_filename", soildepth_input_filename_textbox.Text);
                xwriter.WriteElementString("landuse_input_filename", landuse_input_filename_textbox.Text);
                xwriter.WriteElementString("tillfields_input_filename", tillfields_input_filename_textbox.Text);
                xwriter.WriteElementString("rain_input_filename", rain_input_filename_textbox.Text);
                xwriter.WriteElementString("infil_input_filename", infil_input_filename_textbox.Text);
                xwriter.WriteElementString("evap_input_filename", evap_input_filename_textbox.Text);

                xwriter.WriteElementString("soildepth_constant_value", soildepth_constant_value_box.Text);
                xwriter.WriteElementString("landuse_constant_value", landuse_constant_value_box.Text);
                xwriter.WriteElementString("tillfields_constant_value", tillfields_constant_textbox.Text);
                xwriter.WriteElementString("rain_constant_value", rainfall_constant_value_box.Text);
                xwriter.WriteElementString("infil_constant_value", infil_constant_value_box.Text);
                xwriter.WriteElementString("evap_constant_value", evap_constant_value_box.Text);

                xwriter.WriteElementString("check_fill_sinks_before", XmlConvert.ToString(fill_sinks_before_checkbox.Checked));
                xwriter.WriteElementString("check_fill_sinks_during", XmlConvert.ToString(fill_sinks_during_checkbox.Checked));

                xwriter.WriteEndElement();


                xwriter.WriteStartElement("Run");
                xwriter.WriteElementString("runs_radiobutton", XmlConvert.ToString(runs_checkbox.Checked));
                xwriter.WriteElementString("number_runs", Number_runs_textbox.Text);

                xwriter.WriteStartElement("Specialsettings");
                xwriter.WriteElementString("Spitsbergen", XmlConvert.ToString(Ik_ben_Marijn.Checked));
                xwriter.WriteElementString("Luxembourg", XmlConvert.ToString(version_lux_checkbox.Checked));
                //xwriter.WriteElementString("other", XmlConvert.ToString(runs_checkbox.Checked));
                xwriter.WriteEndElement();

                xwriter.WriteStartElement("CalibrationSensitivity");
                xwriter.WriteElementString("calibration_active_button", XmlConvert.ToString(Calibration_button.Checked));
                xwriter.WriteElementString("calibration_ratios_string", calibration_ratios_textbox.Text);
                xwriter.WriteElementString("calibration_levels", calibration_levels_textbox.Text);
                xwriter.WriteElementString("calibration_ratio_reduction_per_level", calibration_ratio_reduction_parameter_textbox.Text);
                xwriter.WriteEndElement();

                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Output");

                xwriter.WriteStartElement("File_Output");

                xwriter.WriteStartElement("Moment_of_Output");
                xwriter.WriteElementString("final_output_checkbox", XmlConvert.ToString(Final_output_checkbox.Checked));
                xwriter.WriteElementString("regular_output_checkbox", XmlConvert.ToString(Regular_output_checkbox.Checked));
                xwriter.WriteElementString("years_between_outputs", Box_years_output.Text);
                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Type_of_Output");
                xwriter.WriteElementString("cumulative", XmlConvert.ToString(cumulative_output_checkbox.Checked));
                xwriter.WriteElementString("annual", XmlConvert.ToString(annual_output_checkbox.Checked));
                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Maps_to_Output");
                xwriter.WriteElementString("alti", XmlConvert.ToString(Altitude_output_checkbox.Checked));
                xwriter.WriteElementString("altichange", XmlConvert.ToString(Alt_change_output_checkbox.Checked));
                xwriter.WriteElementString("soildepth", XmlConvert.ToString(Soildepth_output_checkbox.Checked));
                xwriter.WriteElementString("all_processes", XmlConvert.ToString(all_process_output_checkbox.Checked));
                xwriter.WriteElementString("waterflow", XmlConvert.ToString(water_output_checkbox.Checked));
                xwriter.WriteElementString("depressions", XmlConvert.ToString(depressions_output_checkbox.Checked));
                xwriter.WriteElementString("diagnostics", XmlConvert.ToString(diagnostic_output_checkbox.Checked));
                xwriter.WriteEndElement();

                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Other_outputs");

                xwriter.WriteStartElement("Avi_Animation");
                xwriter.WriteElementString("Save_as_AVI", XmlConvert.ToString(checkBoxGenerateAVIFile.Checked));
                xwriter.WriteElementString("AVIfilename", textBoxAVIFile.Text);
                xwriter.WriteElementString("years_between_AVI_outputs", saveintervalbox.Text);
                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Google_Animation");
                xwriter.WriteElementString("Save_as_Google", XmlConvert.ToString(googleAnimationCheckbox.Checked));
                xwriter.WriteElementString("Google_filename", googleAnimationTextBox.Text);
                xwriter.WriteElementString("years_between_Google_outputs", googAnimationSaveInterval.Text);
                xwriter.WriteElementString("begindate", googleBeginDate.Text);
                xwriter.WriteElementString("UTM", XmlConvert.ToString(UTMgridcheckbox.Checked));
                xwriter.WriteElementString("South", XmlConvert.ToString(UTMsouthcheck.Checked));
                xwriter.WriteElementString("UTMzone", UTMzonebox.Text);
                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Timeseries");
                xwriter.WriteElementString("total_erosion", XmlConvert.ToString(timeseries.timeseries_total_ero_check.Checked));
                xwriter.WriteElementString("total_deposition", XmlConvert.ToString(timeseries.timeseries_total_dep_check.Checked));
                xwriter.WriteElementString("net_erosion", XmlConvert.ToString(timeseries.timeseries_net_ero_check.Checked));
                xwriter.WriteElementString("SDR", XmlConvert.ToString(timeseries.timeseries_SDR_check.Checked));
                xwriter.WriteElementString("total_average_alt", XmlConvert.ToString(timeseries.timeseries_total_average_alt_check.Checked));
                xwriter.WriteElementString("total_rain", XmlConvert.ToString(timeseries.timeseries_total_rain_check.Checked));
                xwriter.WriteElementString("total_infil", XmlConvert.ToString(timeseries.timeseries_total_infil_check.Checked));
                xwriter.WriteElementString("total_evap", XmlConvert.ToString(timeseries.timeseries_total_evap_check.Checked));
                xwriter.WriteElementString("total_outflow", XmlConvert.ToString(timeseries.timeseries_total_outflow_check.Checked));
                xwriter.WriteElementString("wet_cells", XmlConvert.ToString(timeseries.timeseries_number_waterflow_check.Checked));
                xwriter.WriteElementString("eroded_cells", XmlConvert.ToString(timeseries.timeseries_number_erosion_check.Checked));
                xwriter.WriteElementString("deposited_cells", XmlConvert.ToString(timeseries.timeseries_number_dep_check.Checked));
                xwriter.WriteElementString("cell_altitude", XmlConvert.ToString(timeseries.timeseries_cell_altitude_check.Checked));
                xwriter.WriteElementString("cell_waterflow", XmlConvert.ToString(timeseries.timeseries_cell_waterflow_check.Checked));
                xwriter.WriteElementString("waterflow_threshold", timeseries.timeseries_textbox_waterflow_threshold.Text);
                xwriter.WriteElementString("erosion_threshold", timeseries.timeseries_textbox_erosion_threshold.Text);
                xwriter.WriteElementString("deposition_threshold", timeseries.timeseries_textbox_deposition_threshold.Text);
                xwriter.WriteElementString("cell_row", timeseries.timeseries_textbox_cell_row.Text);
                xwriter.WriteElementString("cell_col", timeseries.timeseries_textbox_cell_col.Text);
                xwriter.WriteElementString("total_OM_input", XmlConvert.ToString(timeseries.total_OM_input_checkbox.Checked));
                xwriter.WriteElementString("total_average_soil_thickness", XmlConvert.ToString(timeseries.total_average_soilthickness_checkbox.Checked));
                xwriter.WriteElementString("total_phys_weath", XmlConvert.ToString(timeseries.total_phys_weath_checkbox.Checked));
                xwriter.WriteElementString("total_chem_weath", XmlConvert.ToString(timeseries.total_chem_weath_checkbox.Checked));
                xwriter.WriteElementString("total_fine_formed", XmlConvert.ToString(timeseries.total_fine_formed_checkbox.Checked));
                xwriter.WriteElementString("total_fine_eluviated", XmlConvert.ToString(timeseries.total_fine_eluviated_checkbox.Checked));
                xwriter.WriteElementString("total_mass_bioturbed", XmlConvert.ToString(timeseries.total_mass_bioturbed_checkbox.Checked));
                xwriter.WriteElementString("timeseries_soil_depth", XmlConvert.ToString(timeseries.timeseries_soil_depth_checkbox.Checked));
                xwriter.WriteElementString("timeseries_soil_mass", XmlConvert.ToString(timeseries.timeseries_soil_mass_checkbox.Checked));
                xwriter.WriteElementString("timeseries_coarser", XmlConvert.ToString(timeseries.timeseries_coarser_checkbox.Checked));
                xwriter.WriteElementString("timeseries_thicker", XmlConvert.ToString(timeseries.timeseries_number_soil_thicker_checkbox.Checked));
                xwriter.WriteElementString("soil_cell", timeseries.timeseries_soil_cell_col.Text);
                xwriter.WriteElementString("soil_col", timeseries.timeseries_soil_cell_row.Text);
                xwriter.WriteElementString("coarser_fraction", timeseries.timeseries_soil_coarser_fraction_textbox.Text);
                xwriter.WriteElementString("thicker_threshold", timeseries.timeseries_soil_thicker_textbox.Text);
                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Profiles");
                xwriter.WriteElementString("profile1_row", XmlConvert.ToString(profile.radio_pro1_row.Checked));
                xwriter.WriteElementString("profile1_col", XmlConvert.ToString(profile.radio_pro1_col.Checked));
                xwriter.WriteElementString("profile2_row", XmlConvert.ToString(profile.radio_pro2_row.Checked));
                xwriter.WriteElementString("profile2_col", XmlConvert.ToString(profile.radio_pro2_col.Checked));
                xwriter.WriteElementString("profile3_row", XmlConvert.ToString(profile.radio_pro3_row.Checked));
                xwriter.WriteElementString("profile3_col", XmlConvert.ToString(profile.radio_pro3_col.Checked));
                xwriter.WriteElementString("p1_number", profile.p1_row_col_box.Text);
                xwriter.WriteElementString("p2_number", profile.p2_row_col_box.Text);
                xwriter.WriteElementString("p3_number", profile.p3_row_col_box.Text);
                xwriter.WriteElementString("p1_waterflow", XmlConvert.ToString(profile.check_waterflow_profile1.Checked));
                xwriter.WriteElementString("p1_altitude", XmlConvert.ToString(profile.check_altitude_profile1.Checked));
                xwriter.WriteElementString("p2_waterflow", XmlConvert.ToString(profile.check_waterflow_profile2.Checked));
                xwriter.WriteElementString("p2_altitude", XmlConvert.ToString(profile.check_altitude_profile2.Checked));
                xwriter.WriteElementString("p3_waterflow", XmlConvert.ToString(profile.check_waterflow_profile3.Checked));
                xwriter.WriteElementString("p3_altitude", XmlConvert.ToString(profile.check_altitude_profile3.Checked));
                xwriter.WriteEndElement();

                xwriter.WriteElementString("maps_required", XmlConvert.ToString(view_maps_checkbox.Checked));

                xwriter.WriteEndElement();

                xwriter.WriteStartElement("Soilfractions");
                xwriter.WriteElementString("coarsefrac", soildata.coarsebox.Text);
                xwriter.WriteElementString("sandfrac", soildata.sandbox.Text);
                xwriter.WriteElementString("siltfrac", soildata.siltbox.Text);
                xwriter.WriteElementString("clayfrac", soildata.claybox.Text);
                xwriter.WriteElementString("fclayfrac", soildata.fineclaybox.Text);
                xwriter.WriteEndElement();
                //End the document
                xwriter.WriteEndDocument();

                //Flush the xml document to the underlying stream and
                //close the underlying stream. The data will not be
                //written out to the stream until either the Flush()
                //method is called or the Close() method is called.
                xwriter.Close();

                this.Text = basetext + " (" + Path.GetFileName(cfgname) + ")";
            }
        }


        #endregion


        #region mapping and drawing code


        void calc_hillshade() // 
        {
            //Local variables
            int row, col;

            double slopemax;
            double slope;
            int slopetot;
            double local_Illumination;

            // Initialize Hillshade Paramaters
            double azimuth = 315 * (3.141592654 / 180); // Default of 315 degrees converted to radians
            double altitude = 45 * (3.141592654 / 180); // Default of 45 degrees converted to radians

            for (row = 1; row < GlobalMethods.nr - 1; row++)
            {
                for (col = 1; col < GlobalMethods.nc - 1; col++)
                {
                    if (GlobalMethods.dtm[row, col] != -9999 && GlobalMethods.dtm[row, col] > 0)
                    {
                        slopemax = 0;
                        slope = 0;
                        slopetot = 0;

                        // Do slope analysis and Aspect Calculation first
                        if (GlobalMethods.dtm[row, col] > GlobalMethods.dtm[row - 1, col] && GlobalMethods.dtm[row - 1, col] != -9999) // North 0
                        {
                            slope = Math.Pow((GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row - 1, col]) / root, 1);
                            if (slope > slopemax)
                            {
                                slopemax = slope;
                                slopetot++;
                                GlobalMethods.aspect[row, col] = 0 * (3.141592654 / 180);
                            }

                        }
                        if (GlobalMethods.dtm[row, col] > GlobalMethods.dtm[row - 1, col + 1] && GlobalMethods.dtm[row - 1, col + 1] != -9999) // Northeast 45
                        {
                            slope = Math.Pow((GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row - 1, col + 1]) / GlobalMethods.d_x, 1);
                            if (slope > slopemax)
                            {
                                slopemax = slope;
                                slopetot++;
                                GlobalMethods.aspect[row, col] = 45 * (3.141592654 / 180);
                            }
                        }
                        if (GlobalMethods.dtm[row, col] > GlobalMethods.dtm[row, col + 1] && GlobalMethods.dtm[row, col + 1] != -9999) // East 90
                        {
                            slope = Math.Pow((GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row, col + 1]) / root, 1);
                            if (slope > slopemax)
                            {
                                slopemax = slope;
                                slopetot++;
                                GlobalMethods.aspect[row, col] = 90 * (3.141592654 / 180);
                            }
                        }
                        if (GlobalMethods.dtm[row, col] > GlobalMethods.dtm[row + 1, col + 1] && GlobalMethods.dtm[row + 1, col + 1] != -9999) // SouthEast 135
                        {
                            slope = Math.Pow((GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + 1, col + 1]) / root, 1);
                            if (slope > slopemax)
                            {
                                slopemax = slope;
                                slopetot++;
                                GlobalMethods.aspect[row, col] = 135 * (3.141592654 / 180);
                            }

                        }
                        if (GlobalMethods.dtm[row, col] > GlobalMethods.dtm[row + 1, col] && GlobalMethods.dtm[row + 1, col] != -9999) // South 180
                        {
                            slope = Math.Pow((GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + 1, col]) / GlobalMethods.d_x, 1);
                            if (slope > slopemax)
                            {
                                slopemax = slope;
                                slopetot++;
                                GlobalMethods.aspect[row, col] = 180 * (3.141592654 / 180);
                            }
                        }
                        if (GlobalMethods.dtm[row, col] > GlobalMethods.dtm[row + 1, col - 1] && GlobalMethods.dtm[row + 1, col - 1] != -9999) // SouthWest 225
                        {
                            slope = Math.Pow((GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row + 1, col - 1]) / root, 1);
                            if (slope > slopemax)
                            {
                                slopemax = slope;
                                slopetot++;
                                GlobalMethods.aspect[row, col] = 225 * (3.141592654 / 180);
                            }
                        }
                        if (GlobalMethods.dtm[row, col] > GlobalMethods.dtm[row, col - 1] && GlobalMethods.dtm[row, col - 1] != -9999) // West 270
                        {
                            slope = Math.Pow((GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row, col - 1]) / root, 1);
                            if (slope > slopemax)
                            {
                                slopemax = slope;
                                slopetot++;
                                GlobalMethods.aspect[row, col] = 270;
                            }
                        }
                        if (GlobalMethods.dtm[row, col] > GlobalMethods.dtm[row - 1, col - 1] && GlobalMethods.dtm[row - 1, col - 1] != -9999) // Northwest 315
                        {
                            slope = Math.Pow((GlobalMethods.dtm[row, col] - GlobalMethods.dtm[row - 1, col - 1]) / GlobalMethods.d_x, 1);
                            if (slope > slopemax)
                            {
                                slopemax = slope;
                                slopetot++;
                                GlobalMethods.aspect[row, col] = 315 * (3.141592654 / 180);
                            }
                        }

                        if (slope > 0) GlobalMethods.slopeAnalysis[row, col] = slopemax;// Tom's: (slope/slopetot); ?

                        // Convert slope to radians
                        GlobalMethods.slopeAnalysis[row, col] = System.Math.Atan(GlobalMethods.slopeAnalysis[row, col]);


                        // Do Hillshade Calculation
                        local_Illumination = 255 * ((System.Math.Cos(azimuth)
                                                     * System.Math.Sin(GlobalMethods.slopeAnalysis[row, col])
                                                     * System.Math.Cos(GlobalMethods.aspect[row, col] - azimuth))
                                                   + (System.Math.Sin(altitude)
                                                     * System.Math.Cos(GlobalMethods.slopeAnalysis[row, col])));

                        GlobalMethods.hillshade[row, col] = System.Math.Abs(local_Illumination);
                    }
                }
            }

        }       // End calc_hillshade() <JOE 20051605- end>

        void Color_HSVtoRGB()   // <JOE 20051605>
        {
            // Convert HSV to RGB.
            // Made this a seperate function as it is called multiple times in draw_map().

            if (sat == 0)
            {
                // If sat is 0, all colors are the same.
                // This is some flavor of gray.
                red = val;
                green = val;
                blue = val;
            }
            else
            {
                double pFactor;
                double qFactor;
                double tFactor;

                double fractionalSector;
                int sectorNumber;
                double sectorPos;

                // The color wheel consists of six 60 degree sectors.
                // Figure out which sector you are in.
                sectorPos = hue / 60;
                sectorNumber = (int)(Math.Floor(sectorPos));

                // get the fractional part of the sector.
                // That is, how many degrees into the sector are you?
                fractionalSector = sectorPos - sectorNumber;

                // Calculate values for the three axes
                // of the color. 
                pFactor = val * (1 - sat);
                qFactor = val * (1 - (sat * fractionalSector));
                tFactor = val * (1 - (sat * (1 - fractionalSector)));

                // Assign the fractional colors to r, g, and b based on the sector the angle is in.
                switch (sectorNumber)
                {
                    case 0:
                        red = val;
                        green = tFactor;
                        blue = pFactor;
                        break;
                    case 1:
                        red = qFactor;
                        green = val;
                        blue = pFactor;
                        break;
                    case 2:
                        red = pFactor;
                        green = val;
                        blue = tFactor;
                        break;
                    case 3:
                        red = pFactor;
                        green = qFactor;
                        blue = val;
                        break;
                    case 4:
                        red = tFactor;
                        green = pFactor;
                        blue = val;
                        break;
                    case 5:
                        red = val;
                        green = pFactor;
                        blue = qFactor;
                        break;
                }
            }
        }

        void draw_map(System.Drawing.Graphics graphics)// <JMW 20041018>
        {
            Debug.WriteLine("\n--drawing maps--");
            Graphics objGraphics;
            objGraphics = Graphics.FromImage(GlobalMethods.m_objDrawingSurface);
            objGraphics.Clear(SystemColors.Control);

            int row, col, z;
            int redcol = 0, greencol = 0, bluecol = 0, alphacol = 255;
            int t = 0;

            // Set Graphics Display Size
            if (GlobalMethods.nc <= 0) GlobalMethods.nc = 1;

            //set scaling of graphics - so X bmp pixels to every model pixel.
            t = graphics_scale;

            // These loop through the entire grid
            // DEM <JOE 20050905>
            if (1 == 1)
            {
                double zDEM;
                double zCalc, zMin = 100000, zMax = 0, zRange, hsMin = 0, hsMax = 255, hs;
                double valMin = 0;
                double valMax = 1;

                calc_hillshade();       // Call up routine 

                // First, find max, min and range of DEM and Hillshade
                for (row = 0; row < GlobalMethods.nr; row++)
                {
                    for (col = 0; col < GlobalMethods.nc; col++)
                    {
                        zCalc = GlobalMethods.dtm[row, col];
                        if (zCalc != -9999)
                        {
                            if (zCalc < zMin) zMin = zCalc;
                            if (zCalc > zMax) zMax = zCalc;
                            hs = GlobalMethods.hillshade[row, col];
                            if (hs < hsMin) hsMin = hs;
                            if (hs > hsMax) hsMax = hs;
                            if (zCalc < -9900) { Debug.WriteLine(" Cell " + row + " " + col + " has altitude " + zCalc); }
                        }



                    }
                }
                if (zMin < 0) { zMin = 0; }
                zRange = zMax - zMin;     //makes the value

                Debug.WriteLine(" zMax, zMin, zRange : " + zMax + " " + zMin + " " + zRange);

                for (row = 0; row < GlobalMethods.nr - 1; row++)
                {
                    for (col = 0; col < GlobalMethods.nc - 1; col++)
                    {
                        if (GlobalMethods.dtm[row, col] > 0)
                        {
                            // HILLSHADE: Draw first underneath
                            // set gray scale intensity
                            hue = 360;  // hue doesn't matter for gray shade
                            sat = 0;        // ensures gray shade
                            valMin = 0;
                            valMax = 1;
                            val = ((GlobalMethods.hillshade[row, col] / 255) * (valMax - valMin)) + valMin; // uses maximum contrast

                            Color_HSVtoRGB();   // Call up color conversion
                            redcol = System.Convert.ToInt32(red * 255);
                            greencol = System.Convert.ToInt32(green * 255);
                            bluecol = System.Convert.ToInt32(blue * 255);
                            alphacol = 255;

                            SolidBrush brush2 = new SolidBrush(Color.FromArgb(alphacol, redcol, greencol, bluecol));
                            objGraphics.FillRectangle(brush2, (col) * t, (row) * t, t, t);

                            // DEM: Colouring based on altitude	
                            zDEM = (GlobalMethods.dtm[row, col]);
                            // Sets hue based on desired color range (in decimal degrees; max 360)
                            double hueMin = 20; //30 = orange   // 20 = reddish orange
                            double hueMax = 140; //85 = green   // 140 = greenish blue
                            hue = hueMax - (((zDEM - zMin) / (zRange)) * (hueMax - hueMin)); // Reverse

                            // Set saturation based on desired range
                            double satMin = 0.50;
                            double satMax = 0.95;
                            sat = (((zDEM - zMin) / (zRange)) * (satMax - satMin)) + satMin;
                            //sat = 0; // Use for grey-scale DEM only!

                            // Set value based on desired range
                            valMin = 0.40;
                            valMax = 0.80;
                            val = (((zDEM - zMin) / (zRange)) * (valMax - valMin)) + valMin;
                            //							val = valMax - (((zDEM - zMin)/(zRange)) * (valMax - valMin));

                            Color_HSVtoRGB();   // Call up color conversion
                            redcol = System.Convert.ToInt32(red * 255);
                            greencol = System.Convert.ToInt32(green * 255);
                            bluecol = System.Convert.ToInt32(blue * 255);
                            alphacol = 125;

                            SolidBrush brush = new SolidBrush(Color.FromArgb(alphacol, redcol, greencol, bluecol));
                            objGraphics.FillRectangle(brush, (col) * t, (row) * t, t, t);
                        }   // Close of Entire Grid Mask
                    }       // Close of Column Loop 
                }           // Close of Row Loop
            }               // Close of DEM check box (menuitem34) 


            // All these loop through just the 'Active Area'
            for (row = 0; row < GlobalMethods.nr - 1; row++)
            {
                for (col = 0; col < GlobalMethods.nc - 1; col++)
                {
                    if (1 > 0) // Index masks out so only 'active cells' shown	was if(GlobalMethods.index[row,col]>-9999)		
                    {
                        // Water Depth
                        if (Menu_map_waterflow.Checked == true && GlobalMethods.waterflow_m3[row, col] > 20)
                        {
                            try { z = (int)(GlobalMethods.waterflow_m3[row, col] * 0.1 * contrastMultiplier); }
                            catch { z = 2147483647; } // this is the maximum integer value.
                            if (z < 0) z = 0;
                            if (z > 255) z = 254;
                            greencol = 255 - z;
                            redcol = z;
                            bluecol = 255;
                            alphacol = 255;
                            SolidBrush brush = new SolidBrush(Color.FromArgb(alphacol, redcol, greencol, bluecol));
                            objGraphics.FillRectangle(brush, (col) * t, (row) * t, t, t);
                        }
                        // Total erosion/ deposition
                        if (Menu_map_total_sediment.Checked == true)
                        {
                            if (GlobalMethods.dtmchange[row, col] < -0.05) //eroding
                            {
                                z = (int)(-GlobalMethods.dtmchange[row, col] * 200 * contrastMultiplier);
                                if (z < 0) z = 0;
                                if (z > 254) z = 254;
                                greencol = 255 - z;
                                redcol = z;
                                SolidBrush brush = new SolidBrush(Color.FromArgb(255, greencol, greencol));
                                objGraphics.FillRectangle(brush, (col) * t, (row) * t, t, t);
                            }
                            if (GlobalMethods.dtmchange[row, col] > 0.01) //depositing
                            {
                                z = (int)(GlobalMethods.dtmchange[row, col] * 1000 * contrastMultiplier);
                                if (z < 0) z = 0;
                                if (z > 254) z = 254;
                                greencol = z;
                                redcol = 255 - z;
                                SolidBrush brush = new SolidBrush(Color.FromArgb(redcol, 255, redcol));
                                objGraphics.FillRectangle(brush, (col) * t, (row) * t, t, t);
                            }
                        }
                        // Water erosion/ deposition
                        if (Menu_map_water_ero.Checked == true)
                        {
                            if (GlobalMethods.sum_water_erosion[row, col] / 25 < -0.05) //eroding
                            {
                                z = (int)(-GlobalMethods.sum_water_erosion[row, col] * 40 * contrastMultiplier);
                                if (z < 0) z = 0;
                                if (z > 254) z = 254;
                                greencol = 255 - z;
                                redcol = z;
                                SolidBrush brush = new SolidBrush(Color.FromArgb(255, greencol, greencol));
                                objGraphics.FillRectangle(brush, (col) * t, (row) * t, t, t);
                            }
                            if (GlobalMethods.sum_water_erosion[row, col] / 25 > 0.05) //depositing
                            {
                                z = (int)(GlobalMethods.sum_water_erosion[row, col] * 40 * contrastMultiplier);
                                if (z < 0) z = 0;
                                if (z > 254) z = 254;
                                greencol = z;
                                redcol = 255 - z;
                                SolidBrush brush = new SolidBrush(Color.FromArgb(redcol, 255, redcol));
                                objGraphics.FillRectangle(brush, (col) * t, (row) * t, t, t);
                            }
                        }
                        // Tillage 
                        if (Menu_map_tillage.Checked == true)
                        {
                            if (GlobalMethods.sum_tillage[row, col] < -0.05) //eroding
                            {
                                z = (int)(-GlobalMethods.sum_tillage[row, col] * 200 * contrastMultiplier);
                                if (z < 0) z = 0;
                                if (z > 254) z = 254;
                                greencol = 255 - z;
                                redcol = z;
                                SolidBrush brush = new SolidBrush(Color.FromArgb(255, greencol, greencol));
                                objGraphics.FillRectangle(brush, (col) * t, (row) * t, t, t);
                            }
                            if (GlobalMethods.sum_tillage[row, col] > 0.05) //depositing
                            {
                                z = (int)(GlobalMethods.sum_tillage[row, col] * 200 * contrastMultiplier);
                                if (z < 0) z = 0;
                                if (z > 254) z = 254;
                                greencol = z;
                                redcol = 255 - z;
                                SolidBrush brush = new SolidBrush(Color.FromArgb(redcol, 255, redcol));
                                objGraphics.FillRectangle(brush, (col) * t, (row) * t, t, t);
                            }
                        }
                        // Creep 
                        if (Menu_map_creep.Checked == true)
                        {
                            if (GlobalMethods.sum_creep_grid[row, col] < -0.05) //eroding
                            {
                                z = (int)(-GlobalMethods.sum_creep_grid[row, col] * 200 * contrastMultiplier);
                                if (z < 0) z = 0;
                                if (z > 254) z = 254;
                                greencol = 255 - z;
                                redcol = z;
                                SolidBrush brush = new SolidBrush(Color.FromArgb(255, greencol, greencol));
                                objGraphics.FillRectangle(brush, (col) * t, (row) * t, t, t);
                            }
                            if (GlobalMethods.sum_creep_grid[row, col] > 0.05) //depositing
                            {
                                z = (int)(GlobalMethods.sum_creep_grid[row, col] * 200 * contrastMultiplier);
                                if (z < 0) z = 0;
                                if (z > 254) z = 254;
                                greencol = z;
                                redcol = 255 - z;
                                SolidBrush brush = new SolidBrush(Color.FromArgb(redcol, 255, redcol));
                                objGraphics.FillRectangle(brush, (col) * t, (row) * t, t, t);
                            }
                        }
                        // Weathering
                        if (Menu_map_weathering.Checked == true && GlobalMethods.sum_biological_weathering[row, col] > 0.001)
                        {
                            z = (int)(GlobalMethods.sum_biological_weathering[row, col] * 1000 * contrastMultiplier);
                            if (z < 0) z = 0;
                            if (z > 255) z = 254;
                            greencol = 255;
                            redcol = z;
                            bluecol = 255 - z;
                            alphacol = 255;
                            SolidBrush brush = new SolidBrush(Color.FromArgb(alphacol, redcol, greencol, bluecol));
                            objGraphics.FillRectangle(brush, (col) * t, (row) * t, t, t);
                        }
                        // Landsliding
                        if (Menu_map_landsliding.Checked == true)
                        {
                            if (GlobalMethods.sum_landsliding[row, col] < -0.1) //eroding
                            {
                                z = (int)(-GlobalMethods.sum_landsliding[row, col] * 40 * contrastMultiplier);
                                if (z < 0) z = 0;
                                if (z > 254) z = 254;
                                greencol = 255 - z;
                                redcol = z;
                                SolidBrush brush = new SolidBrush(Color.FromArgb(255, greencol, greencol));
                                objGraphics.FillRectangle(brush, (col) * t, (row) * t, t, t);
                            }
                            if (GlobalMethods.sum_landsliding[row, col] > 0.1) //depositing
                            {
                                z = (int)(GlobalMethods.sum_landsliding[row, col] * 40 * contrastMultiplier);
                                if (z < 0) z = 0;
                                if (z > 254) z = 254;
                                greencol = z;
                                redcol = 255 - z;
                                SolidBrush brush = new SolidBrush(Color.FromArgb(redcol, 255, redcol));
                                objGraphics.FillRectangle(brush, (col) * t, (row) * t, t, t);
                            }
                        }
                        //Critical rainfall (for landsliding)
                        if (Menu_map_critical_rainfall.Checked == true && GlobalMethods.crrain[row, col] != 99) // no colours for unconditionally stable areas
                        {
                            z = (int)(GlobalMethods.crrain[row, col] * 1000 * contrastMultiplier);
                            if (z < 0) z = 0;
                            if (z > 254) z = 254;
                            greencol = z;
                            redcol = 255 - z;
                            SolidBrush brush = new SolidBrush(Color.FromArgb(redcol, 255, redcol));
                            objGraphics.FillRectangle(brush, (col) * t, (row) * t, t, t);
                        }
                    }           // Close of nodata check for 'active' grid only
                }               // Close of Column Loop
            }                   // Close of Row Loop

            objGraphics.Dispose();
            Mapwindow.Image = GlobalMethods.m_objDrawingSurface;

        }

        #endregion


        #region interface behaviour code

        private void End_button_Click(object sender, EventArgs e)
        {

            this.Close();
        }

        private void UTMgridcheckbox_CheckedChanged(object sender, EventArgs e)
        {
            if (UTMgridcheckbox.Checked) { UTMgroupBox.Visible = true; }
            else { UTMgroupBox.Visible = false; }
        }

        private void Menu_aboutbox_Click(object sender, EventArgs e)
        {
            aboutbox.Visible = true;
        }

        private void timeseries_button_Click(object sender, EventArgs e)
        {
            timeseries.Visible = true;
        }

        private void profiles_button_Click(object sender, EventArgs e)
        {
            profile.Visible = true;
        }

        private void landuse_determinator_button_Click(object sender, EventArgs e)
        {
            landuse_determinator.Visible = true;
        }

        private void Water_ero_checkbox_CheckedChanged(object sender, EventArgs e)
        {
            if (Water_ero_checkbox.Checked == false)
            {
                only_waterflow_checkbox.Enabled = false;
                only_waterflow_checkbox.Checked = false;
            }
            if (Water_ero_checkbox.Checked == true) { only_waterflow_checkbox.Enabled = true; }
        }

        private void Form1_Load(object sender, System.EventArgs e)
        {
            //Mapwindow.Height = this.Height - 180;
            //Mapwindow.Width = this.Width - 200;

            HttpWebRequest req;
            HttpWebResponse res;
            try
            {
                req = (HttpWebRequest)WebRequest.Create("http://www.LORICAmodel.nl/");
                res = (HttpWebResponse)req.GetResponse();
            }
            catch (Exception ex)
            {
                /// do nothing.
            }

            //JMW <20040929 -start>
            this.Text = basetext;
            //DoingGraphics = false;
            //JMW <20040929 - end>

        }

        private void check_cnst_soildepth_CheckedChanged_1(object sender, EventArgs e)
        {
            if (check_space_soildepth.Checked == true) // time can never be true,  because the model calculates soildepth
            {
                soildepth_constant_value_box.Enabled = false;
                soildepth_input_filename_textbox.Enabled = true;
            }
            else
            {
                soildepth_constant_value_box.Enabled = true;
                soildepth_input_filename_textbox.Enabled = false;
            }
        }

        private void check_cnst_landuse_CheckedChanged_1(object sender, EventArgs e)
        {
            if (check_space_landuse.Checked == true)
            {
                landuse_constant_value_box.Enabled = false;
                landuse_input_filename_textbox.Enabled = true;
                check_time_landuse.Checked = false;
            }
            if (check_space_landuse.Checked == false && check_time_landuse.Checked == false)
            {
                landuse_constant_value_box.Enabled = true;
                landuse_input_filename_textbox.Enabled = false;
            }
        }

        private void check_cnst_till_fields_CheckedChanged(object sender, EventArgs e)
        {

            if (check_space_till_fields.Checked == true)
            {
                tillfields_constant_textbox.Enabled = false;
                tillfields_input_filename_textbox.Enabled = true;
                check_time_till_fields.Checked = false;
            }
            if (check_space_till_fields.Checked == false && check_time_till_fields.Checked == false)
            {
                tillfields_constant_textbox.Enabled = true;
                tillfields_input_filename_textbox.Enabled = false;
            }
        }

        private void check_cnst_rain_CheckedChanged_1(object sender, EventArgs e)
        {
            if (check_space_rain.Checked == true)
            {
                rainfall_constant_value_box.Enabled = false;
                rain_input_filename_textbox.Enabled = true;
                check_time_rain.Checked = false;
            }
            if (check_space_rain.Checked == false && check_time_rain.Checked == false)
            {
                rainfall_constant_value_box.Enabled = true;
                rain_input_filename_textbox.Enabled = false;
            }
        }

        private void check_cnst_infil_CheckedChanged(object sender, EventArgs e)
        {
            if (check_space_infil.Checked == true)
            {
                infil_constant_value_box.Enabled = false;
                infil_input_filename_textbox.Enabled = true;
                check_time_infil.Checked = false;
            }
            if (check_space_infil.Checked == false && check_time_infil.Checked == false)
            {
                infil_constant_value_box.Enabled = true;
                infil_input_filename_textbox.Enabled = false;
            }
        }

        private void check_cnst_evap_CheckedChanged(object sender, EventArgs e)
        {
            if (check_space_evap.Checked == true)
            {
                evap_constant_value_box.Enabled = false;
                evap_input_filename_textbox.Enabled = true;
                check_time_evap.Checked = false;
            }
            if (check_space_evap.Checked == false && check_time_evap.Checked == false)
            {
                evap_constant_value_box.Enabled = true;
                evap_input_filename_textbox.Enabled = false;
            }
        }

        private void check_time_landuse_CheckedChanged(object sender, EventArgs e)
        {
            if (check_time_landuse.Checked == true)
            {
                landuse_constant_value_box.Enabled = false;
                landuse_input_filename_textbox.Enabled = true;
                check_space_landuse.Checked = false;
            }
            if (check_space_landuse.Checked == false && check_time_landuse.Checked == false)
            {
                landuse_constant_value_box.Enabled = true;
                landuse_input_filename_textbox.Enabled = false;
            }
        }

        private void check_time_tillage_CheckedChanged(object sender, EventArgs e)
        {
            if (check_time_till_fields.Checked == true) // time can only be true when space is also true
            {

                tillfields_constant_textbox.Enabled = false;
                tillfields_input_filename_textbox.Enabled = true;
                check_space_till_fields.Checked = false;
            }
            if (check_space_till_fields.Checked == false && check_time_till_fields.Checked == false)
            {
                tillfields_constant_textbox.Enabled = true;
                tillfields_input_filename_textbox.Enabled = false;
            }
        }

        private void check_time_rain_CheckedChanged(object sender, EventArgs e)
        {
            if (check_time_rain.Checked == true)
            {
                rainfall_constant_value_box.Enabled = false;
                rain_input_filename_textbox.Enabled = true;
                check_space_rain.Checked = false;
            }
            if (check_space_rain.Checked == false && check_time_rain.Checked == false)
            {
                rainfall_constant_value_box.Enabled = true;
                rain_input_filename_textbox.Enabled = false;
            }
        }

        private void check_time_infil_CheckedChanged(object sender, EventArgs e)
        {
            if (check_time_infil.Checked == true)
            {
                infil_constant_value_box.Enabled = false;
                infil_input_filename_textbox.Enabled = true;
                check_space_infil.Checked = false;
            }
            if (check_space_infil.Checked == false && check_time_infil.Checked == false)
            {
                infil_constant_value_box.Enabled = true;
                infil_input_filename_textbox.Enabled = false;
            }
        }

        private void check_time_evap_CheckedChanged(object sender, EventArgs e)
        {
            if (check_time_evap.Checked == true)
            {
                evap_constant_value_box.Enabled = false;
                evap_input_filename_textbox.Enabled = true;
                check_space_evap.Checked = false;
            }
            if (check_space_evap.Checked == false && check_time_evap.Checked == false)
            {
                evap_constant_value_box.Enabled = true;
                evap_input_filename_textbox.Enabled = false;
            }
        }

        private void button8_Click(object sender, EventArgs e)
        {

            //MessageBox.Show()

            /*"Input filenames are not available when both f(row,col) and f(t) are unchecked. In that case, only single values are input
LORICA will use filename in the following way:

1. When f(row,col) is checked but f(t) is not checked, filename is the ascii grid that will be read.
Example: filename = use.asc ; LORICA will read use.asc

2. When f(row,col) and f(t) are checked, filename is the prefix for 
a series of ascii grid files with the timestep following the prefix. 
Example: filename = use.asc ; LORICA will read use1.asc, use2.asc, use3.asc etc

3. When f(row,col) is not checked, but f(t) is checked, filename is the text file containing 
(spatially uniform) timeseries. The number of values in this file should at least equal 
the number of timesteps in the run. LORICA will start using the first value.
Example: rainfall.asc can look like:
0.67
0.54
0.87
0.70
" */
        }

        private void Menu_map_sediment_Click(object sender, EventArgs e)
        {
            Menu_map_total_sediment.Checked = (!Menu_map_total_sediment.Checked);
            if (Menu_map_total_sediment.Checked)
            {
                Menu_map_waterflow.Checked = false;
                Menu_map_creep.Checked = false;
                Menu_map_tillage.Checked = false;
                Menu_map_water_ero.Checked = false;
                Menu_map_weathering.Checked = false;
                Menu_map_landsliding.Checked = false;
                Menu_map_critical_rainfall.Checked = false;
            }
            updateClick = 1;
        }

        private void Menu_map_waterflow_Click(object sender, EventArgs e)
        {
            if (Water_ero_checkbox.Checked)
            {
                Menu_map_waterflow.Checked = (!Menu_map_waterflow.Checked);
                if (Menu_map_waterflow.Checked)
                {
                    Menu_map_total_sediment.Checked = false;
                    Menu_map_creep.Checked = false;
                    Menu_map_tillage.Checked = false;
                    Menu_map_water_ero.Checked = false;
                    Menu_map_weathering.Checked = false;
                    Menu_map_landsliding.Checked = false;
                    Menu_map_critical_rainfall.Checked = false;
                }
                updateClick = 1;
            }
            else
            {
                MessageBox.Show("water erosion is not active");
            }
        }

        private void Menu_map_tillage_Click(object sender, EventArgs e)
        {
            if (Tillage_checkbox.Checked)
            {
                Menu_map_tillage.Checked = (!Menu_map_tillage.Checked);
                if (Menu_map_tillage.Checked)
                {
                    Menu_map_total_sediment.Checked = false;
                    Menu_map_creep.Checked = false;
                    Menu_map_waterflow.Checked = false;
                    Menu_map_water_ero.Checked = false;
                    Menu_map_weathering.Checked = false;
                    Menu_map_landsliding.Checked = false;
                    Menu_map_critical_rainfall.Checked = false;
                }
                updateClick = 1;
            }
            else
            {
                MessageBox.Show("this process is not active");
            }
        }

        private void Menu_map_water_ero_Click(object sender, EventArgs e)
        {
            if (Water_ero_checkbox.Checked)
            {
                Menu_map_water_ero.Checked = (!Menu_map_water_ero.Checked);
                if (Menu_map_water_ero.Checked)
                {
                    Menu_map_total_sediment.Checked = false;
                    Menu_map_creep.Checked = false;
                    Menu_map_waterflow.Checked = false;
                    Menu_map_tillage.Checked = false;
                    Menu_map_weathering.Checked = false;
                    Menu_map_landsliding.Checked = false;
                    Menu_map_critical_rainfall.Checked = false;
                }
                updateClick = 1;
            }
            else
            {
                MessageBox.Show("this process is not active");
            }
        }

        private void Menu_map_creep_Click(object sender, EventArgs e)
        {
            if (creep_active_checkbox.Checked)
            {
                Menu_map_creep.Checked = (!Menu_map_creep.Checked);
                if (Menu_map_creep.Checked)
                {
                    Menu_map_total_sediment.Checked = false;
                    Menu_map_water_ero.Checked = false;
                    Menu_map_waterflow.Checked = false;
                    Menu_map_tillage.Checked = false;
                    Menu_map_weathering.Checked = false;
                    Menu_map_landsliding.Checked = false;
                    Menu_map_critical_rainfall.Checked = false;
                }
                updateClick = 1;
            }
            else
            {
                MessageBox.Show("this process is not active");
            }
        }

        private void Menu_map_landsliding_Click(object sender, EventArgs e)
        {
            if (Landslide_checkbox.Checked)
            {
                Menu_map_landsliding.Checked = (!Menu_map_landsliding.Checked);
                if (Menu_map_landsliding.Checked)
                {
                    Menu_map_total_sediment.Checked = false;
                    Menu_map_water_ero.Checked = false;
                    Menu_map_waterflow.Checked = false;
                    Menu_map_tillage.Checked = false;
                    Menu_map_creep.Checked = false;
                    Menu_map_weathering.Checked = false;
                    Menu_map_critical_rainfall.Checked = false;
                }
                updateClick = 1;
            }
            else
            {
                MessageBox.Show("this process is not active");
            }
        }

        private void Menu_map_critical_rainfall_Click(object sender, EventArgs e)
        {
            if (Landslide_checkbox.Checked)
            {
                Menu_map_critical_rainfall.Checked = (!Menu_map_critical_rainfall.Checked);
                if (Menu_map_critical_rainfall.Checked)
                {
                    Menu_map_total_sediment.Checked = false;
                    Menu_map_water_ero.Checked = false;
                    Menu_map_waterflow.Checked = false;
                    Menu_map_tillage.Checked = false;
                    Menu_map_creep.Checked = false;
                    Menu_map_weathering.Checked = false;
                    Menu_map_landsliding.Checked = false;
                }
                updateClick = 1;
            }
            else
            {
                MessageBox.Show("this process is not active");
            }
        }

        private void Menu_map_weathering_Click(object sender, EventArgs e)
        {
            if (Biological_weathering_checkbox.Checked)
            {
                Menu_map_weathering.Checked = (!Menu_map_weathering.Checked);
                if (Menu_map_weathering.Checked)
                {
                    Menu_map_total_sediment.Checked = false;
                    Menu_map_water_ero.Checked = false;
                    Menu_map_waterflow.Checked = false;
                    Menu_map_tillage.Checked = false;
                    Menu_map_creep.Checked = false;
                    Menu_map_landsliding.Checked = false;
                    Menu_map_critical_rainfall.Checked = false;
                }
                updateClick = 1;
            }
            else
            {
                MessageBox.Show("this process is not active");
            }
        }

        private void View_tabs_checkbox_CheckedChanged(object sender, EventArgs e)
        {
            if (View_tabs_checkbox.Checked)
            {
                tabControl1.Visible = true;
                graphicToGoogleEarthButton.Visible = false;
                Mapwindow.Visible = false;
            }
            else
            {
                tabControl1.Visible = false;
                graphicToGoogleEarthButton.Visible = true;
                Mapwindow.Visible = true;
            }
        }

        private void dtm_input_filename_textbox_Click(object sender, EventArgs e)
        {
            /*FolderBrowserDialog arnaudsdialog = new FolderBrowserDialog();
            if (arnaudsdialog.ShowDialog() == DialogResult.OK)
            {
                dtm_input_filename_textbox.Text = arnaudsdialog.SelectedPath;
                GlobalMethods.Workdir = arnaudsdialog.SelectedPath;
            } */

            OpenFileDialog openFileDialog1 = new OpenFileDialog();
            openFileDialog1.InitialDirectory = GlobalMethods.Workdir;
            openFileDialog1.Filter = "Ascii grids (*.asc)|*.asc|All files (*.*)|*.*";
            openFileDialog1.FilterIndex = 1;
            openFileDialog1.RestoreDirectory = false;

            dtm_input_filename_textbox.Text = GetDialogFileName(openFileDialog1);
        }

        private void soildepth_input_filename_textbox_TextChanged(object sender, EventArgs e)
        {
            OpenFileDialog openFileDialog1 = new OpenFileDialog();

            openFileDialog1.InitialDirectory = GlobalMethods.Workdir;
            openFileDialog1.Filter = "Ascii grids (*.asc)|*.asc|All files (*.*)|*.*";
            openFileDialog1.FilterIndex = 1;
            openFileDialog1.RestoreDirectory = false;


            soildepth_input_filename_textbox.Text = GetDialogFileName(openFileDialog1);
        }

        private void landuse_input_filename_textbox_TextChanged(object sender, EventArgs e)
        {
            OpenFileDialog openFileDialog1 = new OpenFileDialog();

            openFileDialog1.InitialDirectory = GlobalMethods.Workdir;
            openFileDialog1.Filter = "Ascii grids (*.asc)|*.asc|All files (*.*)|*.*";
            openFileDialog1.FilterIndex = 1;
            openFileDialog1.RestoreDirectory = false;


            landuse_input_filename_textbox.Text = GetDialogFileName(openFileDialog1);
        }

        private void tillfields_input_filename_textbox_TextChanged(object sender, EventArgs e)
        {
            OpenFileDialog openFileDialog1 = new OpenFileDialog();

            openFileDialog1.InitialDirectory = GlobalMethods.Workdir;
            openFileDialog1.Filter = "Ascii grids (*.asc)|*.asc|All files (*.*)|*.*";
            openFileDialog1.FilterIndex = 1;
            openFileDialog1.RestoreDirectory = false;


            tillfields_input_filename_textbox.Text = GetDialogFileName(openFileDialog1);
        }

        private void rain_input_filename_textbox_TextChanged(object sender, EventArgs e)
        {
            OpenFileDialog openFileDialog1 = new OpenFileDialog();

            openFileDialog1.InitialDirectory = GlobalMethods.Workdir;
            openFileDialog1.FilterIndex = 1;
            openFileDialog1.RestoreDirectory = false;

            rain_input_filename_textbox.Text = GetDialogFileName(openFileDialog1);
        }

        private void infil_input_filename_textbox_TextChanged(object sender, EventArgs e)
        {
            OpenFileDialog openFileDialog1 = new OpenFileDialog();

            openFileDialog1.InitialDirectory = GlobalMethods.Workdir;
            openFileDialog1.Filter = "Ascii grids (*.asc)|*.asc|All files (*.*)|*.*";
            openFileDialog1.FilterIndex = 1;
            openFileDialog1.RestoreDirectory = false;

            infil_input_filename_textbox.Text = GetDialogFileName(openFileDialog1);
        }

        private void evap_input_filename_textbox_TextChanged(object sender, EventArgs e)
        {
            OpenFileDialog openFileDialog1 = new OpenFileDialog();

            openFileDialog1.InitialDirectory = GlobalMethods.Workdir;
            openFileDialog1.FilterIndex = 1;
            openFileDialog1.RestoreDirectory = false;

            evap_input_filename_textbox.Text = GetDialogFileName(openFileDialog1);
        }

        private void dailyP_TextChanged(object sender, EventArgs e)
        {
            OpenFileDialog openFileDialog1 = new OpenFileDialog();

            openFileDialog1.InitialDirectory = GlobalMethods.Workdir;
            openFileDialog1.FilterIndex = 1;
            openFileDialog1.RestoreDirectory = false;

            dailyP.Text = GetDialogFileName(openFileDialog1);
        }

        private void dailyET0_TextChanged(object sender, EventArgs e)
        {
            OpenFileDialog openFileDialog1 = new OpenFileDialog();

            openFileDialog1.InitialDirectory = GlobalMethods.Workdir;
            openFileDialog1.FilterIndex = 1;
            openFileDialog1.RestoreDirectory = false;

            dailyET0.Text = GetDialogFileName(openFileDialog1);
        }

        private void dailyD_TextChanged(object sender, EventArgs e)
        {
            OpenFileDialog openFileDialog1 = new OpenFileDialog();

            openFileDialog1.InitialDirectory = GlobalMethods.Workdir;
            openFileDialog1.FilterIndex = 1;
            openFileDialog1.RestoreDirectory = false;

            dailyD.Text = GetDialogFileName(openFileDialog1);
        }

        private void trackBar1_Scroll(object sender, EventArgs e)
        {
            contrastMultiplier = contrastFactor[trackBar1.Value];
            this.view_maps_checkbox.Checked = true;
            draw_map(GlobalMethods.mygraphics);
        }

        private void trackBar2_Scroll(object sender, EventArgs e)
        {
            magnifyValue = zoomFactor[this.trackBar2.Value];
            Mapwindow.setZoom();
        }

        private void comboBox1_SelectedValueChanged(object sender, EventArgs e)
        {
            updateClick = 1;
            this.Refresh();
            this.view_maps_checkbox.Checked = true;
            draw_map(GlobalMethods.mygraphics);
        }

        private void soil_specify_button_Click(object sender, EventArgs e)
        {
            soildata.Visible = true;
        }

        #endregion


        #region top level code



        private void start_run(object sender, System.EventArgs e)
        {
            guiVariables.UpdateFields();
            GlobalMethods.setUp(Ik_ben_Marijn.Checked, 0, 0, guiVariables);
            if (MainSimulation != null)
            {
                throw new MemberAccessException();
            }
            MainSimulation = new Simulation(guiVariables);

            //guiVariables.bulk_density_calc = MainSimulation.bulk_density_calc;


            StartThread = Task.Factory.StartNew(() => {
                MainSimulation.main_loop(sender, e);
            });

            End_button.Enabled = true;
            start_button.Enabled = false;
            
        }

        private void UpdateAllFields()
        {
            this.Invoke(new MethodInvoker(() => {
                this.Refresh(); // tjc to enable graphics to be drawn before sending to AVI

                GlobalMethods.dx = guiVariables.DX;
                rockweath_method.SelectedIndex = guiVariables.Rockweath_method;

                InfoStatusPanel.Text = guiVariables.InfoStatusPanel;
                dtm_input_filename_textbox.Text = guiVariables.DTM_input_filename_textbox;
                Number_runs_textbox.Text = guiVariables.Number_runs_textbox;
                this.ProcessStatusPanel.Text = guiVariables.ProcessStatusPanel;
                calibration_ratios_textbox.Text = guiVariables.Calibration_ratios_textbox;
                googAnimationSaveInterval.Text = guiVariables.GoogAnimationSaveInterval;
                UTMzonebox.Text = guiVariables.UTMzonebox;
                saveintervalbox.Text = guiVariables.Saveintervalbox;
                parameter_m_textbox.Text = guiVariables.Parameter_m_textbox;
                parameter_n_textbox.Text = guiVariables.Parameter_n_textbox;
                parameter_K_textbox.Text = guiVariables.Parameter_K_textbox;
                bio_protection_constant_textbox.Text = guiVariables.Bio_protection_constant_textbox;
                rock_protection_constant_textbox.Text = guiVariables.Rock_protection_constant_textbox;
                selectivity_constant_textbox.Text = guiVariables.Selectivity_constant_textbox;
                erosion_threshold_textbox.Text = guiVariables.Erosion_threshold_textbox;
                parameter_ploughing_depth_textbox.Text = guiVariables.Parameter_ploughing_depth_textbox;
                parameter_tillage_constant_textbox.Text = guiVariables.Parameter_tillage_constant_textbox;
                parameter_diffusivity_textbox.Text = guiVariables.Parameter_diffusivity_textbox;
                parameter_P0_textbox.Text = guiVariables.Parameter_P0_textbox;
                parameter_k1_textbox.Text = guiVariables.Parameter_k1_textbox;
                parameter_k2_textbox.Text = guiVariables.Parameter_k2_textbox;
                parameter_Pa_textbox.Text = guiVariables.Parameter_Pa_textbox;
                Tilting_rate_textbox.Text = guiVariables.Tilting_rate_textbox;
                text_lift_row_less.Text = guiVariables.Text_lift_row_less;
                text_lift_row_more.Text = guiVariables.Text_lift_row_more;
                text_lift_col_less.Text = guiVariables.Text_lift_col_less;
                text_lift_col_more.Text = guiVariables.Text_lift_col_more;
                Uplift_rate_textbox.Text = guiVariables.Uplift_rate_textbox;
                tf_W.Text = guiVariables.Tf_W;
                tf_D.Text = guiVariables.Tf_D;
                tf_growth.Text = guiVariables.Tf_growth;
                tf_age.Text = guiVariables.Tf_age;
                tf_freq.Text = guiVariables.Tf_freq;
                Physical_weath_C1_textbox.Text = guiVariables.Physical_weath_C1_textbox;
                physical_weath_constant1.Text = guiVariables.Physical_weath_constant1;
                physical_weath_constant2.Text = guiVariables.Physical_weath_constant2;
                chem_weath_rate_constant_textbox.Text = guiVariables.Chem_weath_rate_constant_textbox;
                chem_weath_depth_constant_textbox.Text = guiVariables.Chem_weath_depth_constant_textbox;
                chem_weath_specific_coefficient_textbox.Text = guiVariables.Chem_weath_specific_coefficient_textbox;
                specific_area_coarse_textbox.Text = guiVariables.Specific_area_coarse_textbox;
                specific_area_sand_textbox.Text = guiVariables.Specific_area_sand_textbox;
                specific_area_silt_textbox.Text = guiVariables.Specific_area_silt_textbox;
                specific_area_clay_textbox.Text = guiVariables.Specific_area_clay_textbox;
                specific_area_fine_clay_textbox.Text = guiVariables.Specific_area_fine_clay_textbox;
                clay_neoform_constant_textbox.Text = guiVariables.Clay_neoform_constant_textbox;
                clay_neoform_C1_textbox.Text = guiVariables.Clay_neoform_C1_textbox;
                clay_neoform_C2_textbox.Text = guiVariables.Clay_neoform_C2_textbox;
                maximum_eluviation_textbox.Text = guiVariables.Maximum_eluviation_textbox;
                eluviation_coefficient_textbox.Text = guiVariables.Eluviation_coefficient_textbox;
                ct_depth_decay.Text = guiVariables.Ct_depth_decay;
                potential_bioturbation_textbox.Text = guiVariables.Potential_bioturbation_textbox;
                bioturbation_depth_decay_textbox.Text = guiVariables.Bioturbation_depth_decay_textbox;
                carbon_input_textbox.Text = guiVariables.Carbon_input_textbox;
                carbon_depth_decay_textbox.Text = guiVariables.Carbon_depth_decay_textbox;
                carbon_humification_fraction_textbox.Text = guiVariables.Carbon_humification_fraction_textbox;
                carbon_y_decomp_rate_textbox.Text = guiVariables.Carbon_y_decomp_rate_textbox;
                carbon_o_decomp_rate_textbox.Text = guiVariables.Carbon_o_decomp_rate_textbox;
                carbon_y_depth_decay_textbox.Text = guiVariables.Carbon_y_depth_decay_textbox;
                carbon_o_twi_decay_textbox.Text = guiVariables.Carbon_o_twi_decay_textbox;
                carbon_o_depth_decay_textbox.Text = guiVariables.Carbon_o_depth_decay_textbox;
                carbon_y_twi_decay_textbox.Text = guiVariables.Carbon_y_twi_decay_textbox;
                landuse_constant_value_box.Text = guiVariables.Landuse_constant_value_box;
                evap_constant_value_box.Text = guiVariables.Evap_constant_value_box;
                infil_constant_value_box.Text = guiVariables.Infil_constant_value_box;
                rainfall_constant_value_box.Text = guiVariables.Rainfall_constant_value_box;
                temp_constant_value_box.Text = guiVariables.Temp_constant_value_box;
                calibration_ratio_reduction_parameter_textbox.Text = guiVariables.Calibration_ratio_reduction_parameter_textbox;
                soildepth_constant_value_box.Text = guiVariables.Soildepth_constant_value_box;
                Box_years_output.Text = guiVariables.Box_years_output;
                soildepth_input_filename_textbox.Text = guiVariables.Soildepth_input_filename_textbox;
                snow_threshold_textbox.Text = guiVariables.Snow_threshold_textbox;
                snowmelt_factor_textbox.Text = guiVariables.Snowmelt_factor_textbox;
                tillfields_input_filename_textbox.Text = guiVariables.Tillfields_input_filename_textbox;
                landuse_input_filename_textbox.Text= guiVariables.Landuse_input_filename_textbox;
                evap_input_filename_textbox.Text = guiVariables.Evap_input_filename_textbox;
                infil_input_filename_textbox.Text = guiVariables.Infil_input_filename_textbox;
                rain_input_filename_textbox.Text = guiVariables.Rain_input_filename_textbox;
                temp_input_filename_textbox.Text = guiVariables.Temp_input_filename_textbox;
                textBox_ls_trans.Text = guiVariables.TextBox_ls_trans;
                textBox_ls_coh.Text = guiVariables.TextBox_ls_coh;
                textBox_ls_bd.Text = guiVariables.TextBox_ls_bd;
                textBox_ls_ifr.Text = guiVariables.TextBox_ls_ifr;
                total_tillage_statuspanel.Text = guiVariables.Total_tillage_statuspanel;
                dailyP.Text = guiVariables.DailyP;
                dailyD.Text = guiVariables.DailyD;
                dailyT_avg.Text = guiVariables.DailyT_avg;
                dailyT_min.Text = guiVariables.DailyT_min;
                dailyT_max.Text = guiVariables.DailyT_max;
                googleBeginDate.Text = guiVariables.GoogleBeginDate;
                textBoxAVIFile.Text = guiVariables.TextBoxAVIFile;
                ini_CaCO3_content.Text = guiVariables.Ini_CaCO3_content;
                calibration_levels_textbox.Text = guiVariables.Calibration_levels_textbox;
                daily_n.Text = guiVariables.Daily_n;
                latitude_deg.Text = guiVariables.Latitude_deg;
                latitude_min.Text = guiVariables.Latitude_min; 


                Water_ero_checkbox.Checked = guiVariables.Water_ero_checkbox;
                Tillage_checkbox.Checked = guiVariables.Tillage_checkbox;
                Landslide_checkbox.Checked = guiVariables.Landslide_checkbox;
                creep_active_checkbox.Checked = guiVariables.Creep_active_checkbox;
                Biological_weathering_checkbox.Checked = guiVariables.Biological_weathering_checkbox;
                Frost_weathering_checkbox.Checked = guiVariables.Frost_weathering_checkbox;
                tilting_active_checkbox.Checked = guiVariables.Tilting_active_checkbox;
                uplift_active_checkbox.Checked = guiVariables.Uplift_active_checkbox;
                soil_phys_weath_checkbox.Checked = guiVariables.Soil_phys_weath_checkbox;
                soil_chem_weath_checkbox.Checked = guiVariables.Soil_chem_weath_checkbox;
                soil_bioturb_checkbox.Checked = guiVariables.Soil_bioturb_checkbox;
                soil_clay_transloc_checkbox.Checked = guiVariables.Soil_clay_transloc_checkbox;
                soil_carbon_cycle_checkbox.Checked = guiVariables.Soil_carbon_cycle_checkbox;
                UTMsouthcheck.Checked = guiVariables.UTMsouthcheck;
                check_time_till_fields.Checked = guiVariables.Check_time_till_fields;
                check_scaling_daily_weather.Checked = guiVariables.Check_scaling_daily_weather;

                UpdateTimeSeries();
                this.Refresh(); // tjc to enable graphics to be drawn before sending to AVI
            }
            ));
        }
        private void UpdateStatusPannel()
        {
            this.Invoke(new MethodInvoker(() => {
                InfoStatusPanel.Text = guiVariables.InfoStatusPanel;
                out_sed_statuspanel.Text = guiVariables.Out_sed_statuspanel;
            }
            ));
        }
        private void UpdateTimePannel()
        {
            this.Invoke(new MethodInvoker(() => {
                TimeStatusPanel.Text = guiVariables.TimeStatusPanel;
            }
            ));
        }
        private void DelUpdateVariables()
        {
            this.Invoke(new MethodInvoker(() => {
                this.Refresh(); // tjc to enable graphics to be drawn before sending to AVI

                guiVariables.DX = GlobalMethods.dx;
                guiVariables.Rockweath_method = rockweath_method.SelectedIndex;
                try{
                    guiVariables.Mapwindow.Size.Width = this.Mapwindow.Size.Width;
                    guiVariables.Mapwindow.Size.Height = this.Mapwindow.Size.Height;
                    guiVariables.Mapwindow.Location.X = this.Mapwindow.Location.X;
                    guiVariables.Mapwindow.Location.Y = this.Mapwindow.Location.Y;
                }
                catch
                {

                }

                #region TextBoxes
                guiVariables.InfoStatusPanel = InfoStatusPanel.Text;
                guiVariables.DTM_input_filename_textbox = dtm_input_filename_textbox.Text;
                guiVariables.Number_runs_textbox = Number_runs_textbox.Text;
                guiVariables.ProcessStatusPanel = this.ProcessStatusPanel.Text;
                guiVariables.Calibration_ratios_textbox = calibration_ratios_textbox.Text;
                guiVariables.GoogAnimationSaveInterval = googAnimationSaveInterval.Text;
                guiVariables.UTMzonebox = UTMzonebox.Text;
                guiVariables.Saveintervalbox = saveintervalbox.Text;
                guiVariables.Parameter_m_textbox = parameter_m_textbox.Text;
                guiVariables.Parameter_n_textbox = parameter_n_textbox.Text;
                guiVariables.Parameter_K_textbox = parameter_K_textbox.Text;
                guiVariables.Bio_protection_constant_textbox = bio_protection_constant_textbox.Text;
                guiVariables.Rock_protection_constant_textbox = rock_protection_constant_textbox.Text;
                guiVariables.Selectivity_constant_textbox = selectivity_constant_textbox.Text;
                guiVariables.Erosion_threshold_textbox = erosion_threshold_textbox.Text;
                guiVariables.Parameter_ploughing_depth_textbox = parameter_ploughing_depth_textbox.Text;
                guiVariables.Parameter_tillage_constant_textbox = parameter_tillage_constant_textbox.Text;
                guiVariables.Parameter_diffusivity_textbox = parameter_diffusivity_textbox.Text;
                guiVariables.Parameter_P0_textbox = parameter_P0_textbox.Text;
                guiVariables.Parameter_k1_textbox = parameter_k1_textbox.Text;
                guiVariables.Parameter_k2_textbox = parameter_k2_textbox.Text;
                guiVariables.Parameter_Pa_textbox = parameter_Pa_textbox.Text;
                guiVariables.Tilting_rate_textbox = Tilting_rate_textbox.Text;
                guiVariables.Text_lift_row_less = text_lift_row_less.Text;
                guiVariables.Text_lift_row_more = text_lift_row_more.Text;
                guiVariables.Text_lift_col_less = text_lift_col_less.Text;
                guiVariables.Text_lift_col_more = text_lift_col_more.Text;
                guiVariables.Uplift_rate_textbox = Uplift_rate_textbox.Text;
                guiVariables.Tf_W = tf_W.Text;
                guiVariables.Tf_D = tf_D.Text;
                guiVariables.Tf_growth = tf_growth.Text;
                guiVariables.Tf_age = tf_age.Text;
                guiVariables.Tf_freq = tf_freq.Text;
                guiVariables.Physical_weath_C1_textbox = Physical_weath_C1_textbox.Text;
                guiVariables.Physical_weath_constant1 = physical_weath_constant1.Text;
                guiVariables.Physical_weath_constant2 = physical_weath_constant2.Text;
                guiVariables.Chem_weath_rate_constant_textbox = chem_weath_rate_constant_textbox.Text;
                guiVariables.Chem_weath_depth_constant_textbox = chem_weath_depth_constant_textbox.Text;
                guiVariables.Chem_weath_specific_coefficient_textbox = chem_weath_specific_coefficient_textbox.Text;
                guiVariables.Specific_area_coarse_textbox = specific_area_coarse_textbox.Text;
                guiVariables.Specific_area_sand_textbox = specific_area_sand_textbox.Text;
                guiVariables.Specific_area_silt_textbox = specific_area_silt_textbox.Text;
                guiVariables.Specific_area_clay_textbox = specific_area_clay_textbox.Text;
                guiVariables.Specific_area_fine_clay_textbox = specific_area_fine_clay_textbox.Text;
                guiVariables.Clay_neoform_constant_textbox = clay_neoform_constant_textbox.Text;
                guiVariables.Clay_neoform_C1_textbox = clay_neoform_C1_textbox.Text;
                guiVariables.Clay_neoform_C2_textbox = clay_neoform_C2_textbox.Text;
                guiVariables.Maximum_eluviation_textbox = maximum_eluviation_textbox.Text;
                guiVariables.Eluviation_coefficient_textbox = eluviation_coefficient_textbox.Text;
                guiVariables.Ct_depth_decay = ct_depth_decay.Text;
                guiVariables.Potential_bioturbation_textbox = potential_bioturbation_textbox.Text;
                guiVariables.Bioturbation_depth_decay_textbox = bioturbation_depth_decay_textbox.Text;
                guiVariables.Carbon_input_textbox = carbon_input_textbox.Text;
                guiVariables.Carbon_depth_decay_textbox = carbon_depth_decay_textbox.Text;
                guiVariables.Carbon_humification_fraction_textbox = carbon_humification_fraction_textbox.Text;
                guiVariables.Carbon_y_decomp_rate_textbox = carbon_y_decomp_rate_textbox.Text;
                guiVariables.Carbon_o_decomp_rate_textbox = carbon_o_decomp_rate_textbox.Text;
                guiVariables.Carbon_y_depth_decay_textbox = carbon_y_depth_decay_textbox.Text;
                guiVariables.Carbon_o_twi_decay_textbox = carbon_o_twi_decay_textbox.Text;
                guiVariables.Carbon_o_depth_decay_textbox = carbon_o_depth_decay_textbox.Text;
                guiVariables.Carbon_y_twi_decay_textbox = carbon_y_twi_decay_textbox.Text;
                guiVariables.Landuse_constant_value_box = landuse_constant_value_box.Text;
                guiVariables.Evap_constant_value_box = evap_constant_value_box.Text;
                guiVariables.Infil_constant_value_box = infil_constant_value_box.Text;
                guiVariables.Rainfall_constant_value_box = rainfall_constant_value_box.Text;
                guiVariables.Temp_constant_value_box = temp_constant_value_box.Text;
                guiVariables.Calibration_ratio_reduction_parameter_textbox = calibration_ratio_reduction_parameter_textbox.Text;
                guiVariables.Soildepth_constant_value_box = soildepth_constant_value_box.Text;
                guiVariables.Box_years_output = Box_years_output.Text;
                guiVariables.Soildepth_input_filename_textbox = soildepth_input_filename_textbox.Text;
                guiVariables.Snow_threshold_textbox = snow_threshold_textbox.Text;
                guiVariables.Snowmelt_factor_textbox = snowmelt_factor_textbox.Text;
                guiVariables.Tillfields_input_filename_textbox = tillfields_input_filename_textbox.Text;
                guiVariables.Landuse_input_filename_textbox = landuse_input_filename_textbox.Text;
                guiVariables.Evap_input_filename_textbox = evap_input_filename_textbox.Text;
                guiVariables.Infil_input_filename_textbox = infil_input_filename_textbox.Text;
                guiVariables.Rain_input_filename_textbox = rain_input_filename_textbox.Text;
                guiVariables.Temp_input_filename_textbox = temp_input_filename_textbox.Text;
                guiVariables.TextBox_ls_trans = textBox_ls_trans.Text;
                guiVariables.TextBox_ls_coh = textBox_ls_coh.Text;
                guiVariables.TextBox_ls_bd = textBox_ls_bd.Text;
                guiVariables.TextBox_ls_ifr = textBox_ls_ifr.Text;
                guiVariables.Total_tillage_statuspanel = total_tillage_statuspanel.Text;
                guiVariables.DailyP = dailyP.Text;
                guiVariables.DailyD = dailyD.Text;
                guiVariables.DailyT_avg = dailyT_avg.Text;
                guiVariables.DailyT_min = dailyT_min.Text;
                guiVariables.DailyT_max = dailyT_max.Text;
                guiVariables.GoogleBeginDate = googleBeginDate.Text;
                guiVariables.TextBoxAVIFile = textBoxAVIFile.Text;
                guiVariables.Ini_CaCO3_content = ini_CaCO3_content.Text;
                guiVariables.Calibration_levels_textbox = calibration_levels_textbox.Text;
                guiVariables.Daily_n = daily_n.Text;
                guiVariables.Latitude_deg = latitude_deg.Text;
                guiVariables.Latitude_min = latitude_min.Text;
                #endregion

                #region CheckBoxes
                guiVariables.Water_ero_checkbox = Water_ero_checkbox.Checked;
                guiVariables.Tillage_checkbox = Tillage_checkbox.Checked;
                guiVariables.Landslide_checkbox = Landslide_checkbox.Checked;
                guiVariables.Creep_active_checkbox = creep_active_checkbox.Checked;
                guiVariables.Biological_weathering_checkbox = Biological_weathering_checkbox.Checked;
                guiVariables.Frost_weathering_checkbox = Frost_weathering_checkbox.Checked;
                guiVariables.Tilting_active_checkbox = tilting_active_checkbox.Checked;
                guiVariables.Uplift_active_checkbox = uplift_active_checkbox.Checked;
                guiVariables.Soil_phys_weath_checkbox = soil_phys_weath_checkbox.Checked;
                guiVariables.Soil_chem_weath_checkbox = soil_chem_weath_checkbox.Checked;
                guiVariables.Soil_bioturb_checkbox = soil_bioturb_checkbox.Checked;
                guiVariables.Soil_clay_transloc_checkbox = soil_clay_transloc_checkbox.Checked;
                guiVariables.Soil_carbon_cycle_checkbox = soil_carbon_cycle_checkbox.Checked;
                guiVariables.UTMsouthcheck = UTMsouthcheck.Checked;
                guiVariables.Check_time_till_fields = check_time_till_fields.Checked;
                guiVariables.Check_scaling_daily_weather = check_scaling_daily_weather.Checked;
                #endregion

                #region TimeSeries
                guiVariables.Timeseries.Timeseries_soil_thicker_textbox = timeseries.timeseries_soil_thicker_textbox.Text;
                guiVariables.Timeseries.Timeseries_soil_coarser_fraction_textbox = timeseries.timeseries_soil_coarser_fraction_textbox.Text;
                guiVariables.Timeseries.Timeseries_soil_cell_row = timeseries.timeseries_soil_cell_row.Text;
                guiVariables.Timeseries.Timeseries_soil_cell_col = timeseries.timeseries_soil_cell_col.Text;
                guiVariables.Timeseries.Timeseries_erosion_threshold = timeseries.timeseries_erosion_threshold;
                guiVariables.Timeseries.Timeseries_deposition_threshold = timeseries.timeseries_deposition_threshold;
                guiVariables.Timeseries.Timeseries_waterflow_threshold = timeseries.timeseries_waterflow_threshold;
                guiVariables.Timeseries.Timeseries_textbox_cell_row = timeseries.timeseries_textbox_cell_row.Text;
                guiVariables.Timeseries.Timeseries_textbox_cell_col = timeseries.timeseries_textbox_cell_col.Text;

                guiVariables.Timeseries.Total_chem_weath_checkbox = timeseries.total_chem_weath_checkbox.Checked;
                guiVariables.Timeseries.Timeseries_number_soil_thicker_checkbox = timeseries.timeseries_number_soil_thicker_checkbox.Checked;
                guiVariables.Timeseries.Total_average_soilthickness_checkbox = timeseries.total_average_soilthickness_checkbox.Checked;
                guiVariables.Timeseries.Timeseries_soil_depth_checkbox = timeseries.timeseries_soil_depth_checkbox.Checked;
                guiVariables.Timeseries.Timeseries_cell_waterflow_check = timeseries.timeseries_cell_waterflow_check.Checked;
                guiVariables.Timeseries.Timeseries_cell_waterflow_check = timeseries.timeseries_cell_waterflow_check.Checked;
                guiVariables.Timeseries.Total_fine_formed_checkbox = timeseries.total_fine_formed_checkbox.Checked;
                guiVariables.Timeseries.Total_mass_bioturbed_checkbox = timeseries.total_mass_bioturbed_checkbox.Checked;
                guiVariables.Timeseries.Total_OM_input_checkbox = timeseries.total_OM_input_checkbox.Checked;
                guiVariables.Timeseries.Total_fine_eluviated_checkbox = timeseries.total_fine_eluviated_checkbox.Checked;
                guiVariables.Timeseries.Timeseries_net_ero_check = timeseries.timeseries_net_ero_check.Checked;
                guiVariables.Timeseries.Timeseries_number_dep_check = timeseries.timeseries_number_dep_check.Checked;
                guiVariables.Timeseries.Timeseries_cell_altitude_check = timeseries.timeseries_cell_altitude_check.Checked;
                guiVariables.Timeseries.Timeseries_number_erosion_check = timeseries.timeseries_number_erosion_check.Checked;
                guiVariables.Timeseries.Timeseries_number_waterflow_check = timeseries.timeseries_number_waterflow_check.Checked;
                guiVariables.Timeseries.Timeseries_SDR_check = timeseries.timeseries_SDR_check.Checked;
                guiVariables.Timeseries.Timeseries_total_average_alt_check = timeseries.timeseries_total_average_alt_check.Checked;
                guiVariables.Timeseries.Timeseries_total_dep_check = timeseries.timeseries_total_dep_check.Checked;
                guiVariables.Timeseries.Timeseries_total_ero_check = timeseries.timeseries_total_ero_check.Checked;
                guiVariables.Timeseries.Timeseries_total_evap_check = timeseries.timeseries_total_evap_check.Checked;
                guiVariables.Timeseries.Timeseries_total_infil_check = timeseries.timeseries_total_infil_check.Checked;
                guiVariables.Timeseries.Timeseries_total_rain_check = timeseries.timeseries_total_rain_check.Checked;
                guiVariables.Timeseries.Total_phys_weath_checkbox = timeseries.total_phys_weath_checkbox.Checked;
                guiVariables.Timeseries.Timeseries_coarser_checkbox = timeseries.timeseries_coarser_checkbox.Checked;
                guiVariables.Timeseries.Timeseries_soil_mass_checkbox = timeseries.timeseries_soil_mass_checkbox.Checked;
                #endregion

                #region Profile
                guiVariables.Profile.P1_row_col_box = profile.p1_row_col_box.Text;
                guiVariables.Profile.P2_row_col_box = profile.p2_row_col_box.Text;
                guiVariables.Profile.P3_row_col_box = profile.p3_row_col_box.Text;

                guiVariables.Profile.Radio_pro1_col = profile.radio_pro1_col.Checked;
                guiVariables.Profile.Check_altitude_profile1 = profile.check_altitude_profile1.Checked;
                guiVariables.Profile.Check_waterflow_profile1 = profile.check_waterflow_profile1.Checked;
                guiVariables.Profile.Radio_pro1_row = profile.radio_pro1_row.Checked;
                guiVariables.Profile.Radio_pro2_col = profile.radio_pro2_col.Checked;
                guiVariables.Profile.Radio_pro2_row = profile.radio_pro2_row.Checked;
                guiVariables.Profile.Radio_pro3_col = profile.radio_pro3_col.Checked;
                guiVariables.Profile.Radio_pro3_row = profile.radio_pro3_row.Checked;
                #endregion

                #region Landuse Determinator

                guiVariables.Landuse_determinator.LU1_Inf_textbox = landuse_determinator.LU1_Inf_textbox.Text;
                guiVariables.Landuse_determinator.LU1_Evap_textbox = landuse_determinator.LU1_Evap_textbox.Text;
                guiVariables.Landuse_determinator.LU1_Ero_textbox = landuse_determinator.LU1_Ero_textbox.Text;
                guiVariables.Landuse_determinator.LU1_Dep_textbox = landuse_determinator.LU1_Dep_textbox.Text;


                guiVariables.Landuse_determinator.LU2_Inf_textbox = landuse_determinator.LU2_Inf_textbox.Text;
                guiVariables.Landuse_determinator.LU2_Evap_textbox = landuse_determinator.LU2_Evap_textbox.Text;
                guiVariables.Landuse_determinator.LU2_Ero_textbox = landuse_determinator.LU2_Ero_textbox.Text;
                guiVariables.Landuse_determinator.LU2_Dep_textbox = landuse_determinator.LU2_Dep_textbox.Text;

                guiVariables.Landuse_determinator.LU3_Inf_textbox = landuse_determinator.LU3_Inf_textbox.Text;
                guiVariables.Landuse_determinator.LU3_Evap_textbox = landuse_determinator.LU3_Evap_textbox.Text;
                guiVariables.Landuse_determinator.LU3_Ero_textbox = landuse_determinator.LU3_Ero_textbox.Text;
                guiVariables.Landuse_determinator.LU3_Dep_textbox = landuse_determinator.LU3_Dep_textbox.Text;

                guiVariables.Landuse_determinator.LU4_Inf_textbox = landuse_determinator.LU4_Inf_textbox.Text;
                guiVariables.Landuse_determinator.LU4_Evap_textbox = landuse_determinator.LU4_Evap_textbox.Text;
                guiVariables.Landuse_determinator.LU4_Ero_textbox = landuse_determinator.LU4_Ero_textbox.Text;
                guiVariables.Landuse_determinator.LU4_Dep_textbox = landuse_determinator.LU4_Dep_textbox.Text;

                guiVariables.Landuse_determinator.LU5_Inf_textbox = landuse_determinator.LU5_Inf_textbox.Text;
                guiVariables.Landuse_determinator.LU5_Evap_textbox = landuse_determinator.LU5_Evap_textbox.Text;
                guiVariables.Landuse_determinator.LU5_Ero_textbox = landuse_determinator.LU5_Ero_textbox.Text;
                guiVariables.Landuse_determinator.LU5_Dep_textbox = landuse_determinator.LU5_Dep_textbox.Text;

                guiVariables.Landuse_determinator.LU6_Inf_textbox = landuse_determinator.LU6_Inf_textbox.Text;
                guiVariables.Landuse_determinator.LU6_Evap_textbox = landuse_determinator.LU6_Evap_textbox.Text;
                guiVariables.Landuse_determinator.LU6_Ero_textbox = landuse_determinator.LU6_Ero_textbox.Text;
                guiVariables.Landuse_determinator.LU6_Dep_textbox = landuse_determinator.LU6_Dep_textbox.Text;

                guiVariables.Landuse_determinator.LU7_Inf_textbox = landuse_determinator.LU7_Inf_textbox.Text;
                guiVariables.Landuse_determinator.LU7_Evap_textbox = landuse_determinator.LU7_Evap_textbox.Text;
                guiVariables.Landuse_determinator.LU7_Ero_textbox = landuse_determinator.LU7_Ero_textbox.Text;
                guiVariables.Landuse_determinator.LU7_Dep_textbox = landuse_determinator.LU7_Dep_textbox.Text;

                guiVariables.Landuse_determinator.LU8_Inf_textbox = landuse_determinator.LU8_Inf_textbox.Text;
                guiVariables.Landuse_determinator.LU8_Evap_textbox = landuse_determinator.LU8_Evap_textbox.Text;
                guiVariables.Landuse_determinator.LU8_Ero_textbox = landuse_determinator.LU8_Ero_textbox.Text;
                guiVariables.Landuse_determinator.LU8_Dep_textbox = landuse_determinator.LU8_Dep_textbox.Text;

                guiVariables.Landuse_determinator.LU9_Inf_textbox = landuse_determinator.LU9_Inf_textbox.Text;
                guiVariables.Landuse_determinator.LU9_Evap_textbox = landuse_determinator.LU9_Evap_textbox.Text;
                guiVariables.Landuse_determinator.LU9_Ero_textbox = landuse_determinator.LU9_Ero_textbox.Text;
                guiVariables.Landuse_determinator.LU9_Dep_textbox = landuse_determinator.LU9_Dep_textbox.Text;

                guiVariables.Landuse_determinator.LU10_Inf_textbox = landuse_determinator.LU10_Inf_textbox.Text;
                guiVariables.Landuse_determinator.LU10_Evap_textbox = landuse_determinator.LU10_Evap_textbox.Text;
                guiVariables.Landuse_determinator.LU10_Ero_textbox = landuse_determinator.LU10_Ero_textbox.Text;
                guiVariables.Landuse_determinator.LU10_Dep_textbox = landuse_determinator.LU10_Dep_textbox.Text;
                #endregion

                #region SoilData
                guiVariables.Soildata.Coarsebox = soildata.coarsebox.Text;
                guiVariables.Soildata.Sandbox = soildata.sandbox.Text;
                guiVariables.Soildata.Siltbox = soildata.siltbox.Text;
                guiVariables.Soildata.Claybox = soildata.claybox.Text;
                guiVariables.Soildata.Fineclaybox = soildata.fineclaybox.Text;
                #endregion
            }
            ));

        }

        private void UpdateTimeSeries()
        {
            timeseries.timeseries_soil_thicker_textbox.Text = guiVariables.Timeseries.Timeseries_soil_thicker_textbox;
            timeseries.timeseries_soil_coarser_fraction_textbox.Text = guiVariables.Timeseries.Timeseries_soil_coarser_fraction_textbox;
            timeseries.timeseries_soil_cell_row.Text = guiVariables.Timeseries.Timeseries_soil_cell_row;
            timeseries.timeseries_soil_cell_col.Text = guiVariables.Timeseries.Timeseries_soil_cell_col;
            timeseries.timeseries_textbox_cell_row.Text = guiVariables.Timeseries.Timeseries_textbox_cell_row;
            timeseries.timeseries_textbox_cell_col.Text = guiVariables.Timeseries.Timeseries_textbox_cell_col;
            timeseries.timeseries_erosion_threshold = guiVariables.Timeseries.Timeseries_erosion_threshold;
            timeseries.timeseries_deposition_threshold = guiVariables.Timeseries.Timeseries_deposition_threshold;
            timeseries.timeseries_waterflow_threshold = guiVariables.Timeseries.Timeseries_waterflow_threshold;

            timeseries.total_chem_weath_checkbox.Checked = guiVariables.Timeseries.Total_chem_weath_checkbox;
            timeseries.timeseries_number_soil_thicker_checkbox.Checked = guiVariables.Timeseries.Timeseries_number_soil_thicker_checkbox;
            timeseries.total_average_soilthickness_checkbox.Checked = guiVariables.Timeseries.Total_average_soilthickness_checkbox;
            timeseries.timeseries_soil_depth_checkbox.Checked = guiVariables.Timeseries.Timeseries_soil_depth_checkbox;
            timeseries.timeseries_cell_waterflow_check.Checked = guiVariables.Timeseries.Timeseries_cell_waterflow_check;
            timeseries.total_fine_formed_checkbox.Checked = guiVariables.Timeseries.Total_fine_formed_checkbox;
            timeseries.total_mass_bioturbed_checkbox.Checked = guiVariables.Timeseries.Total_mass_bioturbed_checkbox;
            timeseries.total_OM_input_checkbox.Checked = guiVariables.Timeseries.Total_OM_input_checkbox;
            timeseries.total_fine_eluviated_checkbox.Checked = guiVariables.Timeseries.Total_fine_eluviated_checkbox;
            timeseries.timeseries_net_ero_check.Checked = guiVariables.Timeseries.Timeseries_net_ero_check;
            timeseries.timeseries_number_dep_check.Checked = guiVariables.Timeseries.Timeseries_number_dep_check;
            timeseries.timeseries_cell_altitude_check.Checked = guiVariables.Timeseries.Timeseries_cell_altitude_check;
            timeseries.timeseries_number_erosion_check.Checked = guiVariables.Timeseries.Timeseries_number_erosion_check;
            timeseries.timeseries_number_waterflow_check.Checked = guiVariables.Timeseries.Timeseries_number_waterflow_check;
            timeseries.timeseries_SDR_check.Checked = guiVariables.Timeseries.Timeseries_SDR_check;
            timeseries.timeseries_total_average_alt_check.Checked = guiVariables.Timeseries.Timeseries_total_average_alt_check;
            timeseries.timeseries_total_dep_check.Checked = guiVariables.Timeseries.Timeseries_total_dep_check;
            timeseries.timeseries_total_ero_check.Checked = guiVariables.Timeseries.Timeseries_total_ero_check;
            timeseries.timeseries_total_evap_check.Checked = guiVariables.Timeseries.Timeseries_total_evap_check;
            timeseries.timeseries_total_infil_check.Checked = guiVariables.Timeseries.Timeseries_total_infil_check;
            timeseries.timeseries_total_rain_check.Checked = guiVariables.Timeseries.Timeseries_total_rain_check;
            timeseries.total_phys_weath_checkbox.Checked = guiVariables.Timeseries.Total_phys_weath_checkbox;
            timeseries.timeseries_coarser_checkbox.Checked = guiVariables.Timeseries.Timeseries_coarser_checkbox;
            timeseries.timeseries_soil_mass_checkbox.Checked = guiVariables.Timeseries.Timeseries_soil_mass_checkbox;
        }
        private void UpdateProfile()
        {
            profile.p1_row_col_box.Text = guiVariables.Profile.P1_row_col_box;
            profile.p2_row_col_box.Text = guiVariables.Profile.P2_row_col_box;
            profile.p3_row_col_box.Text = guiVariables.Profile.P3_row_col_box;

            profile.radio_pro1_col.Checked = guiVariables.Profile.Radio_pro1_col;
            profile.check_altitude_profile1.Checked = guiVariables.Profile.Check_altitude_profile1;
            profile.check_waterflow_profile1.Checked = guiVariables.Profile.Check_waterflow_profile1;
            profile.radio_pro1_row.Checked = guiVariables.Profile.Radio_pro1_row;
            profile.radio_pro2_col.Checked = guiVariables.Profile.Radio_pro2_col;
            profile.radio_pro2_row.Checked = guiVariables.Profile.Radio_pro2_row;
            profile.radio_pro3_col.Checked = guiVariables.Profile.Radio_pro3_col;
            profile.radio_pro3_row.Checked = guiVariables.Profile.Radio_pro3_row;
        }
        private void UpdateLanduse_determinator()
        {
            landuse_determinator.LU1_Inf_textbox.Text = guiVariables.Landuse_determinator.LU1_Inf_textbox;
            landuse_determinator.LU1_Evap_textbox.Text = guiVariables.Landuse_determinator.LU1_Evap_textbox;
            landuse_determinator.LU1_Ero_textbox.Text = guiVariables.Landuse_determinator.LU1_Ero_textbox;
            landuse_determinator.LU1_Dep_textbox.Text = guiVariables.Landuse_determinator.LU1_Dep_textbox;

            landuse_determinator.LU2_Inf_textbox.Text = guiVariables.Landuse_determinator.LU2_Inf_textbox;
            landuse_determinator.LU2_Evap_textbox.Text = guiVariables.Landuse_determinator.LU2_Evap_textbox;
            landuse_determinator.LU2_Ero_textbox.Text = guiVariables.Landuse_determinator.LU2_Ero_textbox;
            landuse_determinator.LU2_Dep_textbox.Text = guiVariables.Landuse_determinator.LU2_Dep_textbox;

            landuse_determinator.LU3_Inf_textbox.Text = guiVariables.Landuse_determinator.LU3_Inf_textbox;
            landuse_determinator.LU3_Evap_textbox.Text = guiVariables.Landuse_determinator.LU3_Evap_textbox;
            landuse_determinator.LU3_Ero_textbox.Text = guiVariables.Landuse_determinator.LU3_Ero_textbox;
            landuse_determinator.LU3_Dep_textbox.Text = guiVariables.Landuse_determinator.LU3_Dep_textbox;

            landuse_determinator.LU4_Inf_textbox.Text = guiVariables.Landuse_determinator.LU4_Inf_textbox;
            landuse_determinator.LU4_Evap_textbox.Text = guiVariables.Landuse_determinator.LU4_Evap_textbox;
            landuse_determinator.LU4_Ero_textbox.Text = guiVariables.Landuse_determinator.LU4_Ero_textbox;
            landuse_determinator.LU4_Dep_textbox.Text = guiVariables.Landuse_determinator.LU4_Dep_textbox;

            landuse_determinator.LU5_Inf_textbox.Text = guiVariables.Landuse_determinator.LU5_Inf_textbox;
            landuse_determinator.LU5_Evap_textbox.Text = guiVariables.Landuse_determinator.LU5_Evap_textbox;
            landuse_determinator.LU5_Ero_textbox.Text = guiVariables.Landuse_determinator.LU5_Ero_textbox;
            landuse_determinator.LU5_Dep_textbox.Text = guiVariables.Landuse_determinator.LU5_Dep_textbox;

            landuse_determinator.LU6_Inf_textbox.Text = guiVariables.Landuse_determinator.LU6_Inf_textbox;
            landuse_determinator.LU6_Evap_textbox.Text = guiVariables.Landuse_determinator.LU6_Evap_textbox;
            landuse_determinator.LU6_Ero_textbox.Text = guiVariables.Landuse_determinator.LU6_Ero_textbox;
            landuse_determinator.LU6_Dep_textbox.Text = guiVariables.Landuse_determinator.LU6_Dep_textbox;

            landuse_determinator.LU7_Inf_textbox.Text = guiVariables.Landuse_determinator.LU7_Inf_textbox;
            landuse_determinator.LU7_Evap_textbox.Text = guiVariables.Landuse_determinator.LU7_Evap_textbox;
            landuse_determinator.LU7_Ero_textbox.Text = guiVariables.Landuse_determinator.LU7_Ero_textbox;
            landuse_determinator.LU7_Dep_textbox.Text = guiVariables.Landuse_determinator.LU7_Dep_textbox;

            landuse_determinator.LU8_Inf_textbox.Text = guiVariables.Landuse_determinator.LU8_Inf_textbox;
            landuse_determinator.LU8_Evap_textbox.Text = guiVariables.Landuse_determinator.LU8_Evap_textbox;
            landuse_determinator.LU8_Ero_textbox.Text = guiVariables.Landuse_determinator.LU8_Ero_textbox;
            landuse_determinator.LU8_Dep_textbox.Text = guiVariables.Landuse_determinator.LU8_Dep_textbox;

            landuse_determinator.LU9_Inf_textbox.Text = guiVariables.Landuse_determinator.LU9_Inf_textbox;
            landuse_determinator.LU9_Evap_textbox.Text = guiVariables.Landuse_determinator.LU9_Evap_textbox;
            landuse_determinator.LU9_Ero_textbox.Text = guiVariables.Landuse_determinator.LU9_Ero_textbox;
            landuse_determinator.LU9_Dep_textbox.Text = guiVariables.Landuse_determinator.LU9_Dep_textbox;

            landuse_determinator.LU10_Inf_textbox.Text = guiVariables.Landuse_determinator.LU10_Inf_textbox;
            landuse_determinator.LU10_Evap_textbox.Text = guiVariables.Landuse_determinator.LU10_Evap_textbox;
            landuse_determinator.LU10_Ero_textbox.Text = guiVariables.Landuse_determinator.LU10_Ero_textbox;
            landuse_determinator.LU10_Dep_textbox.Text = guiVariables.Landuse_determinator.LU10_Dep_textbox;
        }
        private void UpdateSoildata()
        {
            soildata.coarsebox.Text = guiVariables.Soildata.Coarsebox;
            soildata.sandbox.Text = guiVariables.Soildata.Sandbox;
            soildata.siltbox.Text = guiVariables.Soildata.Siltbox;
            soildata.claybox.Text = guiVariables.Soildata.Claybox;
            soildata.fineclaybox.Text = guiVariables.Soildata.Fineclaybox;
        }
        private string GetDialogFileName()
        {
            OpenFileDialog openFileDialog1 = new OpenFileDialog();

            string FileName = "";
            var t = new Thread((ThreadStart)(() => {
                if (openFileDialog1.ShowDialog() == DialogResult.OK)
                {
                    FileName = openFileDialog1.FileName;
                    //return FileName;
                }
                else return;

            }));

            t.SetApartmentState(ApartmentState.STA);
            t.Start();
            t.Join();

            return FileName;
        }
        private string GetDialogFileName(OpenFileDialog openFileDialog1)
        {
            string FileName = "";

            var t = new Thread((ThreadStart)(() => {
                if (openFileDialog1.ShowDialog() == DialogResult.OK)
                {
                    FileName = openFileDialog1.FileName;
                    //return FileName;
                }
                else return;

            }));

            t.SetApartmentState(ApartmentState.STA);
            t.Start();
            t.Join();

            return FileName;
        }

        private void checkBox1_CheckedChanged(object sender, EventArgs e)
        {

        }

        private void label13_Click(object sender, EventArgs e)
        {

        }

        private void textBox3_TextChanged(object sender, EventArgs e)
        {

        }

        private void statusBar1_PanelClick(object sender, StatusBarPanelClickEventArgs e)
        {

        }

        private void checkBox1_CheckedChanged_1(object sender, EventArgs e)
        {

        }

        private void parameter_conv_textbox_TextChanged(object sender, EventArgs e)
        {

        }

        private void label95_Click_1(object sender, EventArgs e)
        {

        }

        private void label111_Click(object sender, EventArgs e)
        {

        }

        private void label112_Click(object sender, EventArgs e)
        {

        }

        private void textBox3_TextChanged_4(object sender, EventArgs e)
        {

        }

        /*private void Quicksort()        
        {
            //Quicksort(unsorted, 0, unsorted.Length - 1);
            int GlobalMethods.number_of_data_cells;
            int i = 0;
            if (t == 0)  // only in the first timestep;
            {
                GlobalMethods.number_of_data_cells = 0;
                for (row = 0; row < GlobalMethods.nr; row++)  
                {
                    for (col = 0; col < GlobalMethods.nc; col++)
                    {
                        if (GlobalMethods.dtm[row, col] != -9999)
                        {
                            GlobalMethods.index[i] = GlobalMethods.dtm[row, col]; GlobalMethods.row_index[i] = row; GlobalMethods.col_index[i] = col;
                            i++;
                        }
                    }
                }
                GlobalMethods.number_of_data_cells = i;
            }
            i = 0; j = GlobalMethods.number_of_data_cells;           
            IComparable pivot = elements[(left + right) / 2];             
            while (i <= j)            {                
                while (elements[i].CompareTo(pivot) < 0)                
                {                    
                    i++;                
                }                 
                while (elements[j].CompareTo(pivot) > 0)                
                {                    
                    j--;                
                }                 
                if (i <= j)                
                {                    
                    // Swap                    
                    IComparable tmp = elements[i];                    
                    elements[i] = elements[j];                    
                    elements[j] = tmp;                     
                    i++;                    
                    j--;                
                }            
            }             
            // Recursive calls            
            if (left < j)            
            {                
                Quicksort(elements, left, j);            
            }             
            if (i < right)            
            {                
                Quicksort(elements, i, right);            
            }        
        } */

#endregion


    }

    #region for projection of coordinates
    class point
    {
        //lat long variables
        public double xcoord;
        public double ycoord;
        public int UTMzone;
        public bool south;
        double transParallelX = 446.448;
        double transParallelY = -125.157;
        double transParallelZ = 542.060;
        double scaleChange = -20.4894 * 0.000001;
        double rotX = (0.1502 / 3600) * (Math.PI / 180);
        double rotY = (0.2470 / 3600) * (Math.PI / 180);
        double rotZ = (0.8421 / 3600) * (Math.PI / 180);
        double a = 6377563.396; //airy 1830 semi-major axis
        double b = 6356256.910; //airy 1830 semi-minor axis
        double a2 = 6378137.000;
        double b2 = 6356752.3142;
        double eSquared = 0;
        double eSquared2 = 0;
        double nO = -100000;//northing of true origin
        double eO = 400000;//easting of true origin
        double fO = 0.9996012717;//scale factor
        double latTrue = 49.0 * (Math.PI / 180.0);//latitude of true origin
        double longTrue = -2.0 * (Math.PI / 180.0);//longitude of true origin
        double psiHash, MBig, v, v2, v3, nLittle, rho, nSquare = 0;
        double vii, viii, ix, Tx2, xi, xii, xiia = 0;
        double helmertX, helmertY, helmertZ, cartX, cartY, cartZ;
        double Height2 = 0;
        double Height3 = 0;
        double finalLati, finalLongi, latiRad, longiRad = 0;
        double rootXYSqr = 0;
        double PHI1, PHI2, PHI = 0;

        public point(double theXcoord, double theYcoord)//constructor
        {
            this.xcoord = theXcoord;
            this.ycoord = theYcoord;
        }

        public void transformPoint()//british os to lat long
        {
            eSquared = (Math.Pow(a, 2) - Math.Pow(b, 2)) / Math.Pow(a, 2);
            Height2 = 0;
            Height3 = 0;
            //OSGB36 easting and northing to OSGB36 latitude and longitude (lower left corner of DTM)
            psiHash = ((this.ycoord - nO) / (a * fO)) + latTrue;
            nLittle = (a - b) / (a + b);
            MBig = b * fO * (((1 + nLittle + ((5.0 / 4.0) * Math.Pow(nLittle, 2)) + ((5.0 / 4.0) * Math.Pow(nLittle, 3))) * (psiHash -
                     latTrue))
               - (((3.0 * nLittle) + (3.0 * Math.Pow(nLittle, 2)) + ((21.0 / 8.0) * Math.Pow(nLittle, 3))) *
                     Math.Sin(psiHash - latTrue) * Math.Cos(psiHash + latTrue))
               + (((15.0 / 8.0 * Math.Pow(nLittle, 2)) + (15.0 / 8.0 * Math.Pow(nLittle, 3))) * Math.Sin(2.0 * (psiHash - latTrue)) * Math.Cos(2.0 * (psiHash + latTrue)))
               - ((35.0 / 24.0 * Math.Pow(nLittle, 3)) * Math.Sin(3.0 * (psiHash - latTrue)) * Math.Cos(3.0 * (psiHash + latTrue))));
            if (Math.Abs((this.ycoord - nO - MBig)) >= 0.01)
            {
                while (Math.Abs((this.ycoord - nO - MBig)) >= 0.01)
                {
                    psiHash = ((this.ycoord - nO - MBig) / (a * fO)) + psiHash;
                    MBig = b * fO * (((1 + nLittle + ((5.0 / 4.0) * Math.Pow(nLittle, 2)) + ((5.0 / 4.0) * Math.Pow(nLittle, 3))) * (psiHash -
                       latTrue))
                 - (((3.0 * nLittle) + (3.0 * Math.Pow(nLittle, 2)) + ((21.0 / 8.0) * Math.Pow(nLittle, 3))) *
                       Math.Sin(psiHash - latTrue) * Math.Cos(psiHash + latTrue))
                 + (((15.0 / 8.0 * Math.Pow(nLittle, 2)) + (15.0 / 8.0 * Math.Pow(nLittle, 3))) * Math.Sin(2.0 * (psiHash - latTrue)) * Math.Cos(2.0 * (psiHash + latTrue)))
                 - ((35.0 / 24.0 * Math.Pow(nLittle, 3)) * Math.Sin(3.0 * (psiHash - latTrue)) * Math.Cos(3.0 * (psiHash + latTrue))));
                }
            }
            v = a * fO * Math.Pow(1 - eSquared * Math.Pow(Math.Sin(psiHash), 2), -.5);
            rho = a * fO * (1 - eSquared) * Math.Pow(1.0 - eSquared * Math.Pow(Math.Sin(psiHash), 2), -1.5);
            nSquare = v / rho - 1.0;
            vii = (Math.Tan(psiHash)) / (2.0 * rho * v);
            viii = ((Math.Tan(psiHash)) / (24.0 * rho * Math.Pow(v, 3))) * (5 + 3.0 * Math.Pow(Math.Tan(psiHash), 2) + nSquare - 9.0 * (Math.Pow(Math.Tan(psiHash), 2) * nSquare));
            ix = (Math.Tan(psiHash) / ((720.0 * rho * Math.Pow(v, 5)))) * (61 + 90.0 * Math.Pow(Math.Tan(psiHash), 2) + 45.0 * Math.Pow(Math.Tan(psiHash), 4));
            Tx2 = (1.0 / Math.Cos(psiHash)) / v;
            xi = (1.0 / Math.Cos(psiHash)) / (6.0 * Math.Pow(v, 3)) * ((v / rho) + (2.0 * Math.Pow(Math.Tan(psiHash), 2)));
            xii = (1.0 / Math.Cos(psiHash)) / (120.0 * Math.Pow(v, 5)) * (5.0 + (28.0 * Math.Pow(Math.Tan(psiHash), 2)) + (24.0 * Math.Pow(Math.Tan(psiHash), 4)));
            xiia = ((1.0 / Math.Cos(psiHash)) / (5040.0 * Math.Pow(v, 7))) * (61.0 + (662.0 * Math.Pow(Math.Tan(psiHash), 2)) + (1320.0 * Math.Pow(Math.Tan(psiHash), 4)) + (720.0 * Math.Pow(Math.Tan(psiHash), 6)));
            latiRad = psiHash - (vii * Math.Pow((this.xcoord - eO), 2)) + (viii * Math.Pow((this.xcoord - eO), 4)) - (ix * Math.Pow((this.xcoord - eO), 6));
            longiRad = longTrue + (Tx2 * (this.xcoord - eO)) - (xi * (Math.Pow((this.xcoord - eO), 3))) + (xii * (Math.Pow((this.xcoord - eO), 5))) - (xiia * (Math.Pow((this.xcoord - eO), 7)));
            //Debug.WriteLine(latiRad * (180 / Math.PI));
            //Debug.WriteLine(longiRad * (180 / Math.PI));
            //OSGB36 Latitude Longitude Height to OSGB36 Cartesian XYZ
            v2 = a / (Math.Sqrt(1 - (eSquared * ((Math.Pow(Math.Sin(latiRad), 2))))));
            cartX = (v2 + Height2) * Math.Cos(latiRad) * Math.Cos(longiRad);
            cartY = (v2 + Height2) * (Math.Cos(latiRad) * Math.Sin(longiRad));
            cartZ = ((v2 * (1 - eSquared)) + Height2) * Math.Sin(latiRad);
            //Debug.WriteLine();
            //Debug.WriteLine(cartX);
            //Debug.WriteLine(cartY);
            //Helmert Datum Transformation (OSGB36 to WGS84)
            helmertX = cartX + (cartX * scaleChange) - (cartY * rotZ) + (cartZ * rotY) + transParallelX;
            helmertY = (cartX * rotZ) + cartY + (cartY * scaleChange) - (cartZ * rotX) + transParallelY;
            helmertZ = (-1 * cartX * rotY) + (cartY * rotX) + cartZ + (cartZ * scaleChange) + transParallelZ;
            //Debug.WriteLine();
            //Debug.WriteLine(helmertX);
            //Debug.WriteLine(helmertY);
            //WGS84 Cartesian XYZ to WGS84 Latitude, longitude and Ellipsoidal height								
            rootXYSqr = Math.Sqrt((Math.Pow(helmertX, 2)) + (Math.Pow(helmertY, 2)));
            eSquared2 = (Math.Pow(a2, 2) - Math.Pow(b2, 2)) / Math.Pow(a2, 2);
            PHI1 = Math.Atan(helmertZ / (rootXYSqr * (1 - eSquared2)));
            v3 = a2 / (Math.Sqrt(1.0 - (eSquared2 * ((Math.Pow(Math.Sin(PHI1), 2))))));
            PHI2 = Math.Atan((helmertZ + (eSquared2 * v3 * (Math.Sin(PHI1)))) / rootXYSqr);
            while (Math.Abs(PHI1 - PHI2) > 0.000000001)
            {
                PHI1 = PHI2;
                v3 = a2 / (Math.Sqrt(1 - (eSquared2 * ((Math.Pow(Math.Sin(PHI1), 2))))));
                PHI2 = Math.Atan((helmertZ + (eSquared2 * v3 * (Math.Sin(PHI1)))) / rootXYSqr);
            }
            PHI = PHI2;
            finalLati = PHI * (180.0 / Math.PI);
            finalLongi = (Math.Atan(helmertY / helmertX)) * (180.0 / Math.PI);
            this.xcoord = finalLongi;
            this.ycoord = finalLati;
        }

        public void transformUTMPoint()
        {
            //transforms coordinates in UTM WGS84 to lat long
            //requires x, y, zone and north/south
            //the code in this function was found at http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html
            //made by Chuck Taylor
            //tested in xls for points in Poland, Turkey and South Africa

            // The code first calculates TM coordinates from UTM coordinates
            // Then calculates corresponding latitude and longitude in radians
            // Before converting back to degrees

            // first version ArT 12-06-09 

            double footpointlatitude = 0;
            double UTMscalefactor = 0.9996;
            double centralmeridian_deg = 0;
            double centralmeridian_rad = 0;
            double y_ = 0;
            double WGS84_sm_a = 6378137;
            double WGS84_sm_b = 6356752.314;
            double n = (WGS84_sm_a - WGS84_sm_b) / (WGS84_sm_a + WGS84_sm_b);

            this.xcoord = (this.xcoord - 500000) / UTMscalefactor;
            if (this.south) this.ycoord = (this.ycoord - 10000000) / UTMscalefactor;
            else this.ycoord /= UTMscalefactor;

            centralmeridian_deg = -183 + (this.UTMzone * 6);
            centralmeridian_rad = centralmeridian_deg / 180 * Math.PI;

            double alpha = (((WGS84_sm_a + WGS84_sm_b) / 2) * (1 + Math.Pow(n, 2) / 4) + (Math.Pow(n, 4) / 64));
            double beta = (3 * n / 2) + (-27 * Math.Pow(n, 3) / 32) + (269 * Math.Pow(n, 5) / 512);
            double gamma = (21 * Math.Pow(n, 2) / 16) + (-55 * Math.Pow(n, 4) / 32);
            double delta = (151 * Math.Pow(n, 3) / 96) + (-417 * Math.Pow(n, 5) / 128);
            double epsilon = (1097 * Math.Pow(n, 4) / 512);

            y_ = this.ycoord / (alpha);
            footpointlatitude = y_ + (beta * Math.Sin(2 * y_)) + (gamma * Math.Sin(4 * y_)) + (delta * Math.Sin(6 * y_)) + (epsilon * Math.Sin(8 * y_));

            double ep2 = (Math.Pow(WGS84_sm_a, 2) - Math.Pow(WGS84_sm_b, 2)) / Math.Pow(WGS84_sm_b, 2);
            double cf = Math.Cos(footpointlatitude);
            double nuf2 = ep2 * Math.Pow(cf, 2);
            double nf = Math.Pow(WGS84_sm_a, 2) / (WGS84_sm_b * Math.Sqrt(1 + nuf2));

            double tf = Math.Tan(footpointlatitude);
            double tf2 = Math.Pow(tf, 2);
            double tf4 = Math.Pow(tf, 4);

            double x1frac = 1 / (1 * Math.Pow(nf, 1) * cf);
            double x2frac = tf / (2 * Math.Pow(nf, 2));
            double x3frac = 1 / (6 * Math.Pow(nf, 3) * cf);
            double x4frac = tf / (24 * Math.Pow(nf, 4));
            double x5frac = 1 / (120 * Math.Pow(nf, 5) * cf);
            double x6frac = tf / (720 * Math.Pow(nf, 6));
            double x7frac = 1 / (5040 * Math.Pow(nf, 7) * cf);
            double x8frac = tf / (40320 * Math.Pow(nf, 8));

            double x2poly = -1 - nuf2;
            double x3poly = -1 - nuf2 - (2 * tf2);
            double x4poly = 5 + 3 * tf2 + 6 * nuf2 - 6 * tf2 * nuf2 - 3 * Math.Pow(nuf2, 2) - 9 * tf2 * nuf2 * nuf2;
            double x5poly = 5 + 28 * tf2 + 24 * tf4 + 6 * nuf2 + 8 * tf2 * nuf2;
            double x6poly = -61 - 90 * tf2 - 45 * tf4 - 107 * nuf2 + 162 * tf2 * nuf2;
            double x7poly = -61 - 662 * tf2 - 1320 * tf4 - 720 * tf4 * tf2;
            double x8poly = 1385 + 3633 * tf2 + 4095 * tf4 + 1575 * tf4 * tf2;

            double latitude = footpointlatitude
                            + x2frac * x2poly * Math.Pow(this.xcoord, 2)
                            + x4frac * x4poly * Math.Pow(this.xcoord, 4)
                            + x6frac * x6poly * Math.Pow(this.xcoord, 6)
                            + x8frac * x8poly * Math.Pow(this.xcoord, 8);
            double longitude = centralmeridian_rad
                            + x1frac * 1 * Math.Pow(this.xcoord, 1)
                            + x3frac * x3poly * Math.Pow(this.xcoord, 3)
                            + x5frac * x5poly * Math.Pow(this.xcoord, 5)
                            + x7frac * x7poly * Math.Pow(this.xcoord, 7);

            this.ycoord = latitude / Math.PI * 180;
            this.xcoord = longitude / Math.PI * 180;
        }

    }
#endregion


}


