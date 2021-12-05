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
    public static class GlobalMethods
    {
        public const int numberofsinks = 10000;           // run the program once to find out the number of sinks. The exact number and any higher number will do....

        public static int[,]
        drainingoutlet_row = new int[numberofsinks, 5],
        drainingoutlet_col = new int[numberofsinks, 5];


        public static bool input_data_error;
        public static string[] inputheader = new string[6];
        public static int nr,
                            nc,
                            row,
                            col;

        public static int n_texture_classes = 5;
        public static double
                dx, dy,	  		// grid size in both row and col
                xcoord, ycoord,
                d_x;

        public static double[,,] sediment_in_transport_kg,         // sediment mass in kg in transport per texture class
                    litter_kg;                     // Litter contents (Luxembourg case study)

        public static double[,]   // double matrices - these are huge memory-eaters and should be minimized 
                    // they only get that memory later, and only when needed
                    waterflow_m3,        //discharge matrix
                    K_fac,
                    P_fac,
                    dz_ero_m,             //altitude change due to erosion  (negative values)
                    dz_sed_m,             //altitude change due to sedimentation   (positive values)
                    young_SOM_in_transport_kg,
                    old_SOM_in_transport_kg,
                    dz_till_bd,
                    lake_sed_m;         //the thickness of lake sediment

        public static int[,]
                                tillfields,         //fields for tillage 
                                treefall_count,     // count number of tree falls
                                status_map,         //geeft aan of een cel een sink, een zadel, een flat of een top is
                                depression,         //geeft aan of een cel bij een meer hoort, en welk meer
                                landuse,            //landuse in classes
                                slidemap,
                                watsh;              // watershed;

        public static double[,] dtm,                //altitude matrixdouble[,]
                                slopeAnalysis,      //for calculation of hillshade
                                aspect,             //aspect for calculation of hillshade
                                soildepth_m,
                                dtmchange,  	    //change in altitude 
                                dz_soil,
                                evapotranspiration,
                                stslope,		    // matrix with steepest descent local slope [rad]
                                crrain,             // matrix with critical steady state rainfall for landsliding [m/d]
                                camf,               // matrix with number of contributing draining cells, multiple flow [-]
                                T_fac,              // matrix with transmissivity [m/d] values
                                C_fac,              // matrix with combined cohesion [-] values
                                Cs_fac,             // matrix with soil cohesion [kPa] values
                                bulkd,              // matrix with bulk density values [g/cm3]
                                intfr,              // matrix with angle of internal friction values [rad]
                                reserv,
                                ero_slid,
                                cel_dist,
                                sed_slid,
                                sed_bud,
                                dh_slid,
                                bedrock_weathering_m,
                                frost_weathering,
                                sum_biological_weathering,
                                sum_frost_weathering,
                                creep,
                                solif,
                                sum_tillage,
                                sum_water_erosion,
                                sum_landsliding,
                                sum_creep_grid,
                                dz_treefall,        // elevation change by tree fall
                                tpi,            //topographic position GlobalMethods.index
                                sum_solifluction,
                                hornbeam_cover_fraction,   //hornbeam fraction 
                                rain,
                                dtmfill_A,
                                sum_uplift,
                                sum_tilting,
                                till_result,
                                veg,
                                Tau,                //for graphics
                                hillshade,          //for graphics
                                infil,              //infiltration matrix
                                original_dtm;
        public static double[,,]     //3D matrices for properties of soil layers in different x y (x,y,z)
                    layerthickness_m,         // : thickness in m 
                    young_SOM_kg,         // : OM mass in kgrams (per voxel = layer * thickness)
                    old_SOM_kg,         // : OM mass in kgrams (per voxel = layer * thickness) 
                    CO3_kg,   // CaCO3, to track decalcification speed. Does not contribute to texture or soil mass (yet) MM
                    bulkdensity;            // : bulkdensity in kg/m3 (over the voxel = layer * thickness)

        public static bool Ik_ben_Marijn;
        public static int t, t_intervene, number_of_data_cells;

        public static double[] index, depressionlevel = new double[numberofsinks];
        public static string[] rowcol_index;
        public static int[] row_index, col_index;  // for sorting the DEM from high to low
        public static GUIVariables guiVariables;

        public static bool merely_calculating_derivatives;
        public static double[,,,]    //4D matrix for soil texture masses in different x,y and z for t texture classes (x,y,z,t)
                    texture_kg;                //mass in kg (per voxel = layer * thickness)


        public static int[,] OSL_age;            // keeps track of the last moment of surfacing: dim1: [GlobalMethods.nr * GlobalMethods.nc * nlayers * ngrains]; dim2: [5] row, col, layer, depo age, stab age
                                                 // Memory restrictions: length of each dimension cannot exceed 2^31 - 1 units (~2.15*10^9). 
                                                 // Limits on memory size depends on properties of system and settings for simulations

        public static int max_soil_layers = 5;
        public static int ngrains = 100; // For OSL calculations



        static ReaderWriterLock WorkdirRWL = new ReaderWriterLock();
        static string workdir;
        public static string Workdir
        {
            get
            {
                WorkdirRWL.AcquireReaderLock(Timeout.Infinite);
                string temp = workdir;
                WorkdirRWL.ReleaseReaderLock();

                return temp;
            }
            set
            {
                WorkdirRWL.AcquireWriterLock(Timeout.Infinite);
                workdir = value;
                WorkdirRWL.ReleaseWriterLock();
            }
        }


        public static void setUp(bool Ik, int T, int T_intervene, GUIVariables gv)
        {
            Ik_ben_Marijn = Ik;
            t = T;
            t_intervene = T_intervene;
            guiVariables = gv;
        }





        public static void calculate_terrain_derivatives()
        {
            //takes the DTM and calculates key derivatives and writes these to ASCII files
            try { dtm_file(guiVariables.DTM_input_filename_textbox); }
            catch { Debug.WriteLine("could not read DEM for derivative calculation"); }

            //declare rasters and memory
            double[,] ledges, nedges, hedges, hhcliff, hlcliff, slhcliff, sllcliff, terruggedindex, ledgeheight;
            int[,] ledgenames;
            terruggedindex = new double[nr, nc];
            ledges = new double[nr, nc];
            nedges = new double[nr, nc];
            hedges = new double[nr, nc];
            hhcliff = new double[nr, nc];
            hlcliff = new double[nr, nc];
            slhcliff = new double[nr, nc];
            sllcliff = new double[nr, nc];
            ledgeheight = new double[nr, nc];
            ledgenames = new int[nr, nc];
            int runner = 0;

            //Topographic Ruggedness Index (Riley, S.J., DeGloria, S.D., Elliot, R., 1999. A terrain ruggedness index that quantifies topographic heterogeneity. Intermt. J. Sci. 5, 23–27.)
            double sum_squared_difference = 0; int num_nbs = 0;
            try
            {
                for (int row = 0; row < nr; row++)
                {
                    for (int col = 0; col < nc; col++)
                    {
                        if (dtm[row, col] != -9999)
                        {
                            sum_squared_difference = 0;
                            num_nbs = 0;
                            for (int i = (-1); i <= 1; i++)
                            {
                                for (int j = (-1); j <= 1; j++)
                                {
                                    if (!(i == 0 && j == 0) && (row + i) >= 0 && (col + j) >= 0 && (row + i) < nr && (col + j) < nc)
                                    {
                                        if (dtm[row + i, col + j] != -9999)
                                        {
                                            sum_squared_difference += Math.Pow((dtm[row, col] - dtm[row + i, col + j]), 2);
                                            num_nbs++;
                                        }
                                    }
                                }
                            }
                            if (num_nbs == 0) { terruggedindex[row, col] = -9999; }
                            else { terruggedindex[row, col] = Math.Sqrt(sum_squared_difference) * (8 / num_nbs); }
                        }
                    }
                }
                out_double("ruggednessindex.asc", terruggedindex);
                Debug.WriteLine("terrain ruggedness index calculation and storage successfull");
            }
            catch { Debug.WriteLine("terrain ruggedness index calculation or storage failed"); }


            // Properties of possible ledges on the hillslope above and below each cell.
            //We need to ingest ledge positions
            try { read_integer("ledgenames.asc", ledgenames); Debug.WriteLine("ledgenames read successfully"); }
            catch { Debug.WriteLine("ledgenames not found"); }

            //then calculate local properties of the landscape around ledges. We expect that ledge positions may be up to 1 cell wrong.
            double maxcliffheight = 0;
            for (row = 0; row < nr; row++)
            {
                for (col = 0; col < nc; col++)
                {
                    ledges[row, col] = -9999;
                    nedges[row, col] = -9999;
                    hedges[row, col] = -9999;
                    hhcliff[row, col] = -9999;
                    hlcliff[row, col] = -9999;
                    slhcliff[row, col] = -9999;
                    sllcliff[row, col] = -9999;
                    ledgeheight[row, col] = -9999;
                    if (dtm[row, col] != -9999)
                    {
                        ledges[row, col] = 0;
                        nedges[row, col] = 0;
                        hedges[row, col] = 0;
                        hhcliff[row, col] = 0;
                        hlcliff[row, col] = 0;
                        slhcliff[row, col] = 0;
                        sllcliff[row, col] = 0;
                        ledgeheight[row, col] = 0;
                        if (ledgenames[row, col] != -9999)
                        {
                            try
                            {
                                maxcliffheight = 0;
                                for (int i = (-1); i <= 1; i++)
                                {
                                    for (int j = (-1); j <= 1; j++)
                                    {
                                        if (!(i == 0 && j == 0) && row + i >= 0 && col + j >= 0 && row + i < nr && col + j < nc)
                                        {
                                            if (dtm[row + i, col + j] != -9999)
                                            {
                                                for (int ii = (-1); ii <= 1; ii++)
                                                {
                                                    for (int jj = (-1); jj <= 1; jj++)
                                                    {
                                                        if (!(i + ii == 0 && j + jj == 0) && row + i + ii >= 0 && col + j + jj >= 0 && row + i + ii < nr && col + j + jj < nc)
                                                        {
                                                            if (dtm[row + i + ii, col + j + jj] != -9999)
                                                            {
                                                                if (Math.Abs(dtm[row + i + ii, col + j + jj] - dtm[row + i, col + j]) > maxcliffheight) { maxcliffheight = Math.Abs(dtm[row + i + ii, col + j + jj] - dtm[row + i, col + j]); }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                ledgeheight[row, col] = maxcliffheight;
                            }
                            catch { }
                        }
                    }
                }
            }
            Debug.WriteLine("ledgeheights determined");

            for (row = 0; row < nr; row++)
            {
                for (col = 0; col < nc; col++)
                {
                    if (dtm[row, col] != -9999)
                    {
                        if (ledgeheight[row, col] == -9999)
                        {
                            double highestledgeheight = 0;
                            int besti = -4, bestj = -4;
                            int i = -1, j = -1;
                            for (i  = (-1); i <= 1; i++)
                            {
                                for (j  = (-1); j <= 1; j++)
                                {
                                    if (!(i == 0 && j == 0) && row + i >= 0 && col + j >= 0 && row + i < nr && col + j < nc)
                                    {
                                        if (ledgeheight[row + i, col + j] > highestledgeheight)
                                        {
                                            highestledgeheight = ledgeheight[row + i, col + j];
                                        }


                                    }
                                }
                            }
                            ledgeheight[row + i, col + j] = highestledgeheight;
                        }
                    }
                }
            }

            //now, we sort the dtm from high to low and walk through it from high to low to assign ledge properties to 
            comb_sort();
            for (runner = number_of_data_cells - 1; runner >= 0; runner--)
            {     // the index is sorted from low to high values, but flow goes from high to low
                if (index[runner] != -9999)
                {
                    row = row_index[runner]; col = col_index[runner];
                    //Debug.WriteLine("now at row " + row + " col " + col + " alt " + dtm[row, col]);
                    if (ledgenames[row, col] != -9999)
                    {
                        //we are on a ledge. Setting and resetting time
                        hhcliff[row, col] = ledgeheight[row, col];
                        slhcliff[row, col] = hhcliff[row, col] / dx;
                        hedges[row, col]++;
                    }
                    else
                    {
                        double tempslhcliff = 0, steepest = 0, steepness, distance, steepdist = 0;
                        for (int i  = (-1); i <= 1; i++)
                        {
                            for (int j  = (-1); j <= 1; j++)
                            {
                                if (!(i == 0 && j == 0) && (row + i >= 0) && (col + j >= 0) && (row + i < nr) && (col + j < nc))
                                {
                                    if (dtm[row + i, col + j] != -9999)
                                    {
                                        if (dtm[row + i, col + j] > dtm[row, col])
                                        {
                                            if (i == 0 || j == 0) { distance = dx; } else { distance = dx * 1.414; }
                                            steepness = (dtm[row + i, col + j] - dtm[row, col]) / distance;
                                            if (steepness > steepest)
                                            {
                                                //we copy the cliffheight from the steepest neighbour cell
                                                steepdist = distance;
                                                steepest = steepness;
                                                hhcliff[row, col] = hhcliff[row + i, col + j]; tempslhcliff = slhcliff[row + i, col + j];
                                                hedges[row, col] = hedges[row + i, col + j];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        //we now have an updated hhcliff, so we can also update slhcliff
                        slhcliff[row, col] = hhcliff[row, col] / (steepdist + hhcliff[row, col] / tempslhcliff);
                    }
                }
            }
            Debug.WriteLine("downslope variables calculated");

            //now, we walk the other way (from low to high in the DTM). 
            for (runner = 0; runner < number_of_data_cells; runner++)
            {
                if (index[runner] != -9999)
                {
                    row = row_index[runner]; col = col_index[runner];
                    if (ledgenames[row, col] != -9999)
                    {
                        //we are on a ledge. Setting and resetting time
                        hlcliff[row, col] = ledgeheight[row, col];
                        sllcliff[row, col] = hlcliff[row, col] / dx;
                        ledges[row, col]++;
                    }
                    else
                    {
                        double tempsllcliff = 0, steepest = 0, steepness = 0, distance = 0, steepdist = 0;
                        for (int i  = (-1); i <= 1; i++)
                        {
                            for (int j  = (-1); j <= 1; j++)
                            {
                                if (!(i == 0 && j == 0) && row + i >= 0 && col + j >= 0 && row + i < nr && col + j < nc)
                                {
                                    if (dtm[row + i, col + j] != -9999)
                                    {
                                        if (dtm[row + i, col + j] < dtm[row, col])
                                        {
                                            if (i == 0 || j == 0) { distance = dx; } else { distance = dx * 1.414; }
                                            steepness = -(dtm[row + i, col + j] - dtm[row, col]) / distance;
                                            if (steepness > steepest)
                                            {
                                                //we copy the cliffheight from the steepest neighbour cell
                                                steepdist = distance;
                                                steepest = steepness;
                                                hlcliff[row, col] = hlcliff[row + i, col + j]; tempsllcliff = sllcliff[row + i, col + j];
                                                ledges[row, col] = ledges[row + i, col + j];
                                            }
                                        }


                                    }
                                }
                            }
                        }
                        //we now have an updated hhcliff, so we can also update slhcliff
                        sllcliff[row, col] = hlcliff[row, col] / (steepdist + hlcliff[row, col] / tempsllcliff);
                    }
                }
            }
            Debug.WriteLine("upslope variables calculated");

            //finally, add up ledges and hedges to get nedges
            for (row = 0; row < nr; row++)
            {
                for (col = 0; col < nc; col++)
                {
                    if (dtm[row, col] != -9999)
                    {
                        nedges[row, col] = ledges[row, col] + hedges[row, col];
                    }
                }
            }

            //now write all these rasters to ascii:
            out_double("ledgeheights.asc", ledgeheight);
            out_double("nedges.asc", nedges);
            out_double("hedges.asc", hedges);
            out_double("ledges.asc", ledges);
            out_double("hhcliff.asc", hhcliff);
            out_double("hlcliff.asc", hlcliff);
            out_double("slhcliff.asc", slhcliff);
            out_double("sllcliff.asc", sllcliff);
            Debug.WriteLine("variables exported to ASCII");
        }

        public static void dtm_file(string name1)
        {

            string FILE_NAME = name1;
            int z, dem_integer_error = 1;
            string[] lineArray2;
            int sp;
            //MessageBox.Show("Opening DEM" + FILE_NAME);
            //MessageBox.Show("Directory " + Directory.GetCurrentDirectory() );

            if (!File.Exists(FILE_NAME))
            {
                MessageBox.Show("No such DEM data file..");
                input_data_error = true;
                return;
            }

            try
            {

                //read headers
                StreamReader sr = File.OpenText(FILE_NAME);
                for (z = 1; z <= 6; z++)
                {
                    inputheader[z - 1] = sr.ReadLine();
                    //MessageBox.Show(inputheader[z - 1]);
                }
                sr.Close();

                // get nc, nr and dx from input headers

                lineArray2 = inputheader[0].Split(new char[] { ' ' });
                sp = 1;
                while (lineArray2[sp] == "") sp++;
                nc = int.Parse(lineArray2[sp]);

                lineArray2 = inputheader[1].Split(new char[] { ' ' });
                sp = 1;
                while (lineArray2[sp] == "") sp++;
                nr = int.Parse(lineArray2[sp]);

                lineArray2 = inputheader[2].Split(new char[] { ' ' });
                sp = 1;
                while (lineArray2[sp] == "") sp++;
                xcoord = double.Parse(lineArray2[sp]);

                lineArray2 = inputheader[3].Split(new char[] { ' ' });
                sp = 1;
                while (lineArray2[sp] == "") sp++;
                ycoord = double.Parse(lineArray2[sp]);

                lineArray2 = inputheader[4].Split(new char[] { ' ' });
                sp = 1;
                while (lineArray2[sp] == "") sp++;
                dx = double.Parse(lineArray2[sp]);

                Debug.WriteLine("read DEM: nr = " + nr + " nc = " + nc);
            }
            catch (Exception ex)
            {
                MessageBox.Show("There is a problem with the header of the DEM file");
                input_data_error = true;
                return;

            }
            int ok = makematrices();
            if (ok == 1)
            { // we have now succesfully made memory reservations for all data layers in the model 

            }
            else
            {
                MessageBox.Show("There is not enough memory for LORICA to run with these settings");
            }
            {

                int col, row, colcounter;
                String input;
                double tttt = 0.00;

                // load dem

                if (!File.Exists(FILE_NAME))
                {
                    MessageBox.Show("No such DEM data file..");
                    input_data_error = true;
                    return;
                }

                StreamReader sr = File.OpenText(FILE_NAME);

                //now skip over the headers.
                for (z = 1; z <= 6; z++)
                {
                    input = sr.ReadLine();
                }
                row = 0;
                while ((input = sr.ReadLine()) != null)  // so not until nr is reached, but until the file is empty
                {
                    string[] lineArray;
                    lineArray = input.Split(new char[] { ' ' });   // so we split the string that we read (readline) from file into an array of strings that each contain a number
                    col = 0;
                    for (colcounter = 0; colcounter <= (lineArray.Length - 1); colcounter++)  // the length of LineArray should equal nc, and therefore run from 0 to nc-1
                    {

                        if (lineArray[colcounter] != "" && col < nc) // but just to make sure, col counts only the non-empty strings in LineArrary (handy for instance when files are double-spaced)
                        {
                            tttt = double.Parse(lineArray[colcounter]);
                            if (Ik_ben_Marijn) { original_dtm[row, col] = tttt; dtm[row, col] = -9999; }
                            else { dtm[row, col] = tttt; }
                            col++;
                            if (double.Parse(lineArray[colcounter]) - Math.Round(double.Parse(lineArray[colcounter])) != 0)
                            {
                                dem_integer_error = 0;
                            }
                        }
                    }
                    row++;


                }
                sr.Close();
                if (dem_integer_error == 1) { MessageBox.Show("Warning: Digital Elevation Model may only contain integer values\n LORICA can proceed, but may experience problems"); }

            }
        }

        public static void out_double(string name4, double[,] output)
        {
            int nn, row, col;
            string FILENAME = name4;
            using (StreamWriter sw = new StreamWriter(FILENAME))
            {
                sw.Write("ncols         " + nc);
                sw.Write("\r\n");
                sw.Write("nrows         " + nr);
                sw.Write("\r\n");
                for (nn = 2; nn <= 5; nn++)
                {
                    sw.Write(inputheader[nn]); sw.Write("\r\n");
                    //MessageBox.Show(inputheader[nn]);
                }
                for (row = 0; row < nr; row++)
                {
                    for (col = 0; col < nc; col++)
                    {
                        sw.Write("{0:F6}", output[row, col]);
                        sw.Write(" ");

                    }
                    sw.Write("\r\n");
                }
                sw.Close();
            }

        } //end out_double

        public static void comb_sort()      //sorts the data cells in a dtm in order of increasing altitude
        {
            // comb sorting by Wlodek Dobosiewicz in 1980
            // http://en.wikipedia.org/wiki/Comb_sort
            // LORICA adaptation by Arnaud Temme june 2009
            //Debug.WriteLine("sorting. nr " + nr + " nc " + nc + " t " + t);
            guiVariables.InfoStatusPanel = "sorting";
            int i = 0;
            double dtm_temp = 0;
            int row_temp = 0, col_temp = 0;
            string rowcol_temp;
            //Debug.WriteLine("sorting. nr " + nr + " nc " + nc + " t " + t);
            if (t == t_intervene)  // only in the first timestep;
            {
                //Debug.WriteLine("normal sorting. nr " + nr + " nc " + nc + " t " + t);
                number_of_data_cells = 0;
                for (row = 0; row < nr; row++)  // why not do this only in the first timestep? And use the existing one as input in subsequent timesteps?
                {
                    for (col = 0; col < nc; col++)
                    {
                        if (dtm[row, col] != -9999)
                        {
                            index[i] = dtm[row, col];
                            row_index[i] = row;
                            col_index[i] = col;
                            rowcol_index[i] = row.ToString() + "." + col.ToString();
                            i++;
                        }
                    }
                }
                number_of_data_cells = i;
            }
            else
            {
                //Debug.WriteLine("alternative sorting. nr " + nr + " nc " + nc + " t " + t);
                for (i = 0; i < number_of_data_cells; i++)
                {
                    index[i] = dtm[row_index[i], col_index[i]];     //merely update the existing index with the adapted altitudes and then sort     
                }
            }
            //displayonscreen(0, 0);
            guiVariables.InfoStatusPanel = "data cells: " + number_of_data_cells;
            //Debug.WriteLine("\n--sorting overview--");
            //Debug.WriteLine("Sorting " + number_of_data_cells + " cells");
            long gap = number_of_data_cells;
            bool swaps;
            long total_swaps = 0;
            //while (gap > 1 && swaps == true)  // in freak? situations, swaps may be false for gap = x, but true for subsequent values of gap
            while (gap > 1)
            {
                if (gap > 1)
                {
                    if (gap == 2) { gap = 1; }
                    gap = Convert.ToInt64(gap / 1.2);
                }
                i = 0;
                swaps = false;
                //this.guiVariables.InfoStatusPanel = "i " + i + " gap " + gap + " tot swaps " + total_swaps;
                //Debug.WriteLine("i " + i + " gap " + gap + " tot swaps " + total_swaps);
                while (i + gap < number_of_data_cells)
                {
                    //if (gap == Convert.ToInt64(number_of_data_cells / 1.2) && i < 10) {Debug.WriteLine("    i " + i + " gap " + gap + " tot swaps " + total_swaps + " alt1 " + index[i] + " (" + row_index[i] + "," + col_index[i] + ") alt2 " + index[i+gap] + " (" + row_index[i+gap] + "," + col_index[i+gap] + ")"); }
                    if (index[i] > index[i + gap])
                    {
                        dtm_temp = index[i]; index[i] = index[i + gap]; index[i + gap] = dtm_temp;
                        row_temp = row_index[i]; row_index[i] = row_index[i + gap]; row_index[i + gap] = row_temp;
                        col_temp = col_index[i]; col_index[i] = col_index[i + gap]; col_index[i + gap] = col_temp;
                        rowcol_temp = rowcol_index[i]; rowcol_index[i] = rowcol_index[i + gap]; rowcol_index[i + gap] = rowcol_temp;
                        swaps = true;
                        total_swaps++;
                    } // end if
                    i++;
                }  // end while
                   //if (gap < 4) { Debug.WriteLine("i " + i + " gap " + gap + " tot swaps " + total_swaps); }
            } //end while
            int sorting_error = 0;
            for (i = 0; i < number_of_data_cells - 1; i++)
            {
                if (index[i] > index[i + 1]) { sorting_error = 1; }
            }
            if (sorting_error == 1)
            {
                Debug.WriteLine(" Sorting error in comb_sort ");
            }
            else
            {
                //Debug.WriteLine(" Sorting test successful ");
            }
        }


        public static int makematrices()
        {
            // Debug.WriteLine("assigning memory");
            // status grids
            if (Ik_ben_Marijn == true) { original_dtm = new double[nr, nc]; }
            dtm = new double[nr, nc];
            if (merely_calculating_derivatives == false)
            {
                OSL_age = new int[nr * nc * max_soil_layers * ngrains, 5];
                soildepth_m = new double[nr, nc];
                dtmchange = new double[nr, nc];
                dz_soil = new double[nr, nc];
                // climate grids
                if (guiVariables.Check_space_evap) { evapotranspiration = new double[nr, nc]; }
                if (guiVariables.Check_space_infil) { infil = new double[nr, nc]; }
                if (guiVariables.Check_space_rain) { rain = new double[nr, nc]; }
                veg = new double[nr, nc];
                // categorical grids
                if (guiVariables.Check_space_landuse) { landuse = new int[nr, nc]; }
            }
            status_map = new int[nr, nc];
            //sorting arrays
            index = new double[nr * nc];
            row_index = new int[nr * nc];
            col_index = new int[nr * nc];
            rowcol_index = new string[nr * nc];
            //others
            depression = new int[nr, nc];
            dtmfill_A = new double[nr, nc];
            if (merely_calculating_derivatives == false)
            {
                if (1 == 1)
                {
                    texture_kg = new double[nr, nc, max_soil_layers, n_texture_classes];    //: mass in kg (per voxel = layer * thickness)
                    layerthickness_m = new double[nr, nc, max_soil_layers];        // : thickness in m 
                    young_SOM_kg = new double[nr, nc, max_soil_layers];         // : OM mass in kg (per voxel = layer * thickness)
                    old_SOM_kg = new double[nr, nc, max_soil_layers];
                    bulkdensity = new double[nr, nc, max_soil_layers];            // : bulkdensity in kg/m3 (over the voxel = layer * thickness)
                }

                if (guiVariables.Water_ero_checkbox)
                {
                    //doubles
                    waterflow_m3 = new double[nr, nc];
                    if (!guiVariables.Only_waterflow_checkbox)
                    {
                        K_fac = new double[nr, nc];
                        P_fac = new double[nr, nc];
                        sediment_in_transport_kg = new double[nr, nc, n_texture_classes];
                        young_SOM_in_transport_kg = new double[nr, nc];
                        old_SOM_in_transport_kg = new double[nr, nc];
                        sum_water_erosion = new double[nr, nc];
                        dz_ero_m = new double[nr, nc];
                        dz_sed_m = new double[nr, nc];
                        lake_sed_m = new double[nr, nc];
                        //depressionsum_texture_kg = new double[n_texture_classes];

                    }

                }
                if (guiVariables.Tillage_checkbox)
                {
                    till_result = new double[nr, nc];
                    sum_tillage = new double[nr, nc];
                    tillfields = new int[nr, nc];
                    dz_till_bd = new double[nr, nc];
                }

                if (guiVariables.Treefall_checkbox)
                {
                    treefall_count = new int[nr, nc];
                    dz_treefall = new double[nr, nc];
                }

                if (guiVariables.Version_lux_checkbox)
                {
                    tpi = new double[nr, nc];
                    hornbeam_cover_fraction = new double[nr, nc];
                    litter_kg = new double[nr, nc, 2];
                }

                if (guiVariables.Solifluction_checkbox)
                {
                    solif = new double[nr, nc];
                    sum_solifluction = new double[nr, nc];
                }
                if (guiVariables.Creep_active_checkbox)
                {
                    creep = new double[nr, nc];
                    sum_creep_grid = new double[nr, nc];
                }
                if (guiVariables.Landslide_checkbox)
                {
                    //doubles
                    stslope = new double[nr, nc];
                    crrain = new double[nr, nc];
                    camf = new double[nr, nc];
                    T_fac = new double[nr, nc];
                    C_fac = new double[nr, nc];
                    Cs_fac = new double[nr, nc];
                    bulkd = new double[nr, nc];
                    intfr = new double[nr, nc];
                    reserv = new double[nr, nc];
                    ero_slid = new double[nr, nc];
                    cel_dist = new double[nr, nc];
                    sed_slid = new double[nr, nc];
                    sed_bud = new double[nr, nc];
                    dh_slid = new double[nr, nc];
                    sum_landsliding = new double[nr, nc];
                    //integers
                    slidemap = new int[nr, nc];
                    watsh = new int[nr, nc];
                }
                if (guiVariables.Biological_weathering_checkbox)
                {
                    bedrock_weathering_m = new double[nr, nc];
                    sum_biological_weathering = new double[nr, nc];
                }
                if (guiVariables.Frost_weathering_checkbox)
                {
                    frost_weathering = new double[nr, nc];
                    sum_frost_weathering = new double[nr, nc];
                }
                if (guiVariables.Tilting_active_checkbox)
                {
                    sum_tilting = new double[nr, nc];
                }
                if (guiVariables.Uplift_active_checkbox)
                {
                    sum_uplift = new double[nr, nc];
                }
                if (guiVariables.Decalcification_checkbox)
                {
                    CO3_kg = new double[nr, nc, max_soil_layers];
                }

            }
            aspect = new double[nr, nc];
            slopeAnalysis = new double[nr, nc];
            hillshade = new double[nr, nc];
            Tau = new double[nr, nc];
            // Debug.WriteLine("memory assigned succesfully");
            return 1;
        }

        public static void read_double(string name2, double[,] map1)
        {
            string FILE_NAME = name2;
            string input;
            double tttt = 0.00;
            int x, y, xcounter;
            if (!File.Exists(FILE_NAME))
            {
                MessageBox.Show("No such double data file " + FILE_NAME);
                GlobalMethods.input_data_error = true;
                return;
            }

            StreamReader sr = File.OpenText(FILE_NAME);

            //read headers
            for (int z = 1; z <= 6; z++)
            {
                input = sr.ReadLine();
            }
            y = 0;

            while ((input = sr.ReadLine()) != null)
            {
                string[] lineArray;
                lineArray = input.Split(new char[] { ' ' });
                xcounter = 0;
                for (x = 0; x <= (lineArray.Length - 1); x++)
                {

                    if (lineArray[x] != "" && xcounter < nc)
                    {


                        try
                        {
                            tttt = double.Parse(lineArray[x]);
                        }
                        catch
                        {
                            MessageBox.Show("Incorrect content " + lineArray[x] + " in file " + FILE_NAME);
                            GlobalMethods.input_data_error = true;
                            return;
                        }
                        map1[y, xcounter] = tttt;
                        xcounter++;
                    }
                }
                y++;

            }
            sr.Close();
        } // end read_double()

        public static void read_integer(string name2, int[,] map1)
        {
            string FILE_NAME = name2;
            string input;
            int tttt = 0;
            int x, y, xcounter;
            Debug.WriteLine(" Reading " + FILE_NAME + " from " + Directory.GetCurrentDirectory());
            if (!File.Exists(FILE_NAME))
            {
                MessageBox.Show("No such data file " + FILE_NAME);
                GlobalMethods.input_data_error = true;
                return;
            }
            StreamReader sr = File.OpenText(FILE_NAME);

            //read headers
            for (int z = 1; z <= 6; z++)
            {
                input = sr.ReadLine();
                /*if (z == 1)
                {
                    string[] lineArray;
                    lineArray = input.Split(new char[] { ' ' });
                    Debug.WriteLine(input + " here " + lineArray[1] + " there " );
                    if (int.Parse(lineArray[1]) != nc)
                    {
                        Debug.WriteLine(filename + " has different cols than the DEM ");
                    }
                }
                if (z == 2)
                {
                    string[] lineArray;
                    lineArray = input.Split(new char[] { ' ' });
                    Debug.WriteLine(lineArray[1]);
                    if (int.Parse(lineArray[1]) != nr)
                    {
                        Debug.WriteLine(filename + " has different rows than the DEM ");
                    }
                } */
            }
            y = 0;
            while ((input = sr.ReadLine()) != null)
            {
                string[] lineArray;
                lineArray = input.Split(new char[] { ' ' });
                xcounter = 0;
                for (x = 0; x <= (lineArray.Length - 1); x++)
                {

                    if (lineArray[x] != "" && xcounter < nc)
                    {
                        try
                        {
                            tttt = int.Parse(lineArray[x]);
                        }
                        catch
                        {
                            MessageBox.Show("Incorrect content " + lineArray[x] + " in file " + FILE_NAME);
                            GlobalMethods.input_data_error = true;
                            return;
                        }
                        map1[y, xcounter] = tttt;
                        xcounter++;
                    }
                }
                y++;

            }
            sr.Close();
            //Debug.WriteLine("completed reading file" + FILE_NAME);
        } // end read_integer()

        public static void read_record(string filename, int[] record)
        {
            string FILE_NAME = filename;
            string input;
            int tttt = 0;
            int y;
            if (!File.Exists(FILE_NAME))
            {
                MessageBox.Show("No such data file " + FILE_NAME);
                GlobalMethods.input_data_error = true;
                return;
            }
            // Debug.WriteLine("reading " + filename + " into record ");
            StreamReader sr = File.OpenText(FILE_NAME);

            //read first line: number of timesteps
            input = sr.ReadLine();
            y = 0;
            int recordsize = 0;
            try { recordsize = System.Convert.ToInt32(input); }
            catch
            {
                MessageBox.Show("Wrong value " + input + " in first line of record " + FILE_NAME);
                GlobalMethods.input_data_error = true;
                return;
            }


            // Debug.WriteLine("reading " + filename + " into record of size " + record.Length);

            // the record size is read from the first line and not necessarily equal to the number of timesteps. 
            // Runs will start from beginning of record and repeat when necessary
            while ((input = sr.ReadLine()) != null)
            {
                if (y >= recordsize)
                {
                    MessageBox.Show("record " + FILE_NAME + " contains more values than expected. Extras are ignored");
                    break;
                }

                try { tttt = int.Parse(input); }
                catch
                {
                    MessageBox.Show("Incorrect content " + input + " in file " + FILE_NAME);
                    GlobalMethods.input_data_error = true;
                    return;
                }
                record[y] = tttt;
                //Debug.WriteLine("value " + y + " in record is " + record[y]);
                y++;
            }
            sr.Close();

        }

        public static void out_mf(string name4, double[,,] output)
        {
            int nn, row, col;
            string FILENAME = name4;
            using (StreamWriter sw = new StreamWriter(FILENAME))
            {
                sw.Write("In n1 n2 n3 n4 n5 n6 n7 n8");
                sw.Write("\r\n");
                for (row = 0; row < nr; row++)
                {
                    for (col = 0; col < nc; col++)
                    {
                        for (int dir = 0; dir < 9; dir++)
                        {
                            sw.Write(guiVariables.OFy_m[row, col, dir]);
                            sw.Write(" ");
                        }
                        sw.Write("\r\n");
                    }
                }

                sw.Write("ncols         " + nc);
                sw.Write("\r\n");
                sw.Write("nrows         " + nr);
                sw.Write("\r\n");

                sw.Close();
            }
        }

        public static void out_integer(string name4, int[,] output)
        {
            int nn, row, col;
            string FILENAME = name4;
            using (StreamWriter sw = new StreamWriter(FILENAME))
            {
                for (nn = 0; nn <= 5; nn++)
                {
                    sw.Write(GlobalMethods.inputheader[nn]); sw.Write("\n");
                }
                for (row = 0; row < nr; row++)
                {
                    for (col = 0; col < nc; col++)
                    {

                        sw.Write(output[row, col]);
                        sw.Write(" ");
                    }

                    sw.Write("\n");
                }
                sw.Close();
            }
        } //end out_integer

        public static void out_profile(string name5, double[,] output, bool row_is_fixed, int row_or_col)
        {
            // WVG 20-10-2010 output a profile file for benefit glorious model of LORICA
            int row, col;
            string FILENAME = name5;
            using (StreamWriter sw = new StreamWriter(FILENAME))

                try
                {
                    if (row_is_fixed)
                    {
                        try
                        {
                            for (col = 0; col < nc; col++)// WVG the number of columns is equal to nc
                            {
                                sw.Write(output[row_or_col, col]);
                                sw.Write(" ");
                            }
                            sw.Write("\n");
                        }
                        catch { Debug.WriteLine("out_profile: error "); }
                    }
                    else  // apparently column is fixed
                    {
                        try
                        {
                            for (row = 0; row < nr; row++)// WVG the number of columns is equal to nc
                            {
                                sw.Write(output[row, row_or_col]);
                                sw.Write(" ");
                            }
                            sw.Write("\n");
                        }
                        catch { Debug.WriteLine("out_profile: error "); }
                    }
                    sw.Close();
                }
                catch { Debug.WriteLine("Profile could not be written"); }

        } //WVG end out_profile

        public static void writesoil(int row, int col)
        {
            int layer;
            double cumthick, midthick;
            string FILENAME = string.Format("{0}\\t{1}_r{2}_c{3}_out_soil.csv", Workdir, t + 1, row, col);
            using (StreamWriter sw = new StreamWriter(FILENAME))
            {
                sw.Write("row, col, t, cumth_m, thick_m, midthick_m, coarse_kg, sand_kg, silt_kg, clay_kg, fine_kg, YOM_kg, OOM_kg, YOM/OOM, f_coarse, f_sand, f_silt, f_clay, f_fineclay");
                sw.Write("\r\n");
                cumthick = 0;
                midthick = 0;
                int t_out = t + 1;
                for (layer = 0; layer < max_soil_layers; layer++) // only the top layer
                {
                    if (layerthickness_m[row, col, layer] > 0)
                    {
                        cumthick += layerthickness_m[row, col, layer];
                        midthick += layerthickness_m[row, col, layer] / 2;
                        double totalweight = texture_kg[row, col, layer, 0] + texture_kg[row, col, layer, 1] + texture_kg[row, col, layer, 2] + texture_kg[row, col, layer, 3] + texture_kg[row, col, layer, 4] + young_SOM_kg[row, col, layer] + old_SOM_kg[row, col, layer];
                        sw.Write(row + "," + col + "," + t_out + "," + cumthick + "," + layerthickness_m[row, col, layer] + "," + midthick + "," + texture_kg[row, col, layer, 0] + "," + texture_kg[row, col, layer, 1] + "," + texture_kg[row, col, layer, 2] + "," + texture_kg[row, col, layer, 3] + "," + texture_kg[row, col, layer, 4] + "," + young_SOM_kg[row, col, layer] + "," + old_SOM_kg[row, col, layer] + "," + young_SOM_kg[row, col, layer] / old_SOM_kg[row, col, layer] + "," + texture_kg[row, col, layer, 0] / totalweight + "," + texture_kg[row, col, layer, 1] / totalweight + "," + texture_kg[row, col, layer, 2] / totalweight + "," + texture_kg[row, col, layer, 3] / totalweight + "," + texture_kg[row, col, layer, 4] / totalweight);
                        sw.Write("\r\n");
                        midthick += layerthickness_m[row, col, layer] / 2;
                    }

                }
                sw.Close();
            }
        }// end writesoil

        public static void writeallsoils()
        {
            int layer;
            double cumthick, midthick, z_layer;
            string FILENAME = string.Format("{0}\\t{1}_out_allsoils.csv", Workdir, t + 1);
            using (StreamWriter sw = new StreamWriter(FILENAME))
            {
                sw.Write("row, col, t, nlayer, cumth_m, thick_m, midthick_m, z, coarse_kg, sand_kg, silt_kg, clay_kg, fine_kg, YOM_kg, OOM_kg, YOM/OOM, f_coarse, f_sand, f_silt, f_clay, f_fineclay, ftotal_clay, f_OM, BD");
                sw.Write("\r\n");
                int t_out = t + 1;
                for (int row = 0; row < nr; row++)
                {
                    for (int col = 0; col < nc; col++)
                    {
                        if (dtm[row, col] != -9999)
                        {
                            cumthick = 0;
                            midthick = 0;
                            z_layer = dtm[row, col];
                            for (layer = 0; layer < max_soil_layers; layer++) // only the top layer
                            {
                                if (layerthickness_m[row, col, layer] > 0)
                                {
                                    cumthick += layerthickness_m[row, col, layer];
                                    midthick += layerthickness_m[row, col, layer] / 2;
                                    double totalweight = texture_kg[row, col, layer, 0] + texture_kg[row, col, layer, 1] + texture_kg[row, col, layer, 2] + texture_kg[row, col, layer, 3] + texture_kg[row, col, layer, 4] + young_SOM_kg[row, col, layer] + old_SOM_kg[row, col, layer];
                                    double totalweight_tex = texture_kg[row, col, layer, 0] + texture_kg[row, col, layer, 1] + texture_kg[row, col, layer, 2] + texture_kg[row, col, layer, 3] + texture_kg[row, col, layer, 4];
                                    sw.Write(row + "," + col + "," + t_out + "," + layer + "," + cumthick + "," + layerthickness_m[row, col, layer] + "," + midthick + "," + z_layer + "," + texture_kg[row, col, layer, 0] + "," + texture_kg[row, col, layer, 1] + "," + texture_kg[row, col, layer, 2] + "," + texture_kg[row, col, layer, 3] + "," + texture_kg[row, col, layer, 4] + "," + young_SOM_kg[row, col, layer] + "," + old_SOM_kg[row, col, layer] + "," + young_SOM_kg[row, col, layer] / old_SOM_kg[row, col, layer] + "," + texture_kg[row, col, layer, 0] / totalweight_tex + "," + texture_kg[row, col, layer, 1] / totalweight_tex + "," + texture_kg[row, col, layer, 2] / totalweight_tex + "," + texture_kg[row, col, layer, 3] / totalweight_tex + "," + texture_kg[row, col, layer, 4] / totalweight_tex + "," + (texture_kg[row, col, layer, 3] + texture_kg[row, col, layer, 4]) / totalweight_tex + "," + (young_SOM_kg[row, col, layer] + old_SOM_kg[row, col, layer]) / (young_SOM_kg[row, col, layer] + old_SOM_kg[row, col, layer] + totalweight_tex) + "," + bulkdensity[row, col, layer]);
                                    sw.Write("\r\n");
                                    midthick += layerthickness_m[row, col, layer] / 2;
                                    z_layer -= layerthickness_m[row, col, layer];
                                }

                            }

                        }
                    }
                }
                sw.Close();
            }

        }// end writeallsoils

        public static void write_longitudinal_profile(int startrow, int startcol, string name4)
        {

            // writes a longitudinal steepest-descent profile starting from a given begin r,c 
            double altidiff, non_lake_altidiff, maxaltidiff;
            int row, col, i, j, non_lake_maxi, non_lake_maxj, maxi, maxj, profilesize = 1000, step;
            double[] profile;
            profile = new double[1000];
            row = startrow;
            col = startcol;
            step = 0;
            string FILENAME = name4;

            while (row - 1 >= 0 && row + 1 < nr && col - 1 >= 0 && col + 1 < nc)
            {   // as long as we have not reached the edge
                Debug.WriteLine("profile now at row %d col %d, alt %.4f\n", row, col, dtm[row, col]);
                altidiff = -9999; non_lake_altidiff = 0; maxaltidiff = 0; non_lake_maxi = 0; non_lake_maxj = 0; maxi = 0; maxj = 0;
                for (i = -1; i <= 1; i++)
                {
                    for (j = -1; j <= 1; j++)
                    {
                        altidiff = dtm[row, col] - dtm[row + i, col + j];
                        //Debug.WriteLine(" profile : nb %d %d, diff %.3f\n",row+i,col+j,altidiff);
                        if (altidiff > maxaltidiff)
                        {
                            maxaltidiff = altidiff; maxi = i; maxj = j;
                        }
                        if (altidiff > non_lake_altidiff && depression[row + i, col + j] == 0)
                        {
                            non_lake_altidiff = altidiff; non_lake_maxi = i; non_lake_maxj = j;
                        }
                    }
                }
                Debug.WriteLine("profile : found lowest nb at %d %d, diff %.3f\n", row + maxi, col + maxj, maxaltidiff);
                if (non_lake_altidiff != 0) { row += non_lake_maxi; col += non_lake_maxj; } //avoid depressions if you can and prevent from falling back
                else { row += maxi; col += maxj; }
                if (maxi == 0 && maxj == 0 || maxaltidiff == 0)
                {
                    Debug.WriteLine("warning : no profile-progress due to sink?\n"); // go straight to the outlet of this depression and count the number of cells in between
                    //break;
                    maxi = drainingoutlet_row[depression[row, col], 0] - row;
                    maxj = drainingoutlet_col[depression[row, col], 0] - col;
                    for (i = 1; i <= (Math.Abs(maxi) + Math.Abs(maxj)); i++)
                    {
                        profile[step] = dtm[row, col] + (i / (Math.Abs(maxi) + Math.Abs(maxj))) * (depressionlevel[depression[row, col]] - dtm[row, col]);
                        Debug.WriteLine("profile %d now %.6f\n", step, profile[step]);
                        step++;
                    }
                    row += maxi;
                    col += maxj;
                }
                else
                {
                    profile[step] = dtm[row, col];
                    step++;
                    if (step > profilesize - 3) { Debug.WriteLine("warning : profilerecord may be too small\n"); break; }
                }
            }

            using (StreamWriter sw = new StreamWriter(FILENAME))
            {
                for (i = 0; i < step + 1; i++)
                {
                    Debug.WriteLine(profile[i]);
                    sw.Write("{0:F6}", profile[i]);
                }
                sw.Write("\r\n");
            }
            Debug.WriteLine("profile contains %d values\n", step);
        } // end write_profile() 

        public static void write_full_output(string filecore, int rows, int cols, int layers, int t)
        {
            int layer, row, col;
            string filename = Workdir + "\\" + filecore + t + ".lrc";
            //Debug.WriteLine("attempting to write output " + filename + " at t " + t);
            using (StreamWriter sw = new StreamWriter(filename))
            {
                try
                {
                    sw.WriteLine("Lorica output header");
                    sw.WriteLine("year " + t);
                    sw.WriteLine("years " + guiVariables.End_time + " every " + int.Parse(guiVariables.Box_years_output));
                    sw.WriteLine("rows " + rows + " cellsize " + dx + " yllcoord " + ycoord);
                    sw.WriteLine("cols " + cols + " xllcoord " + xcoord);
                    sw.WriteLine("layers " + layers);
                    sw.WriteLine("properties 10");
                    sw.WriteLine("propnames elevation thickness_m density_kg_m3 coarse_kg sand_kg silt_kg clay_kg fineclay_kg youngom_kg oldom_kg");
                    sw.WriteLine("Lorica output content");
                }
                catch { Debug.WriteLine(" issue with writing the header of the full output file for this timestep"); }
                try
                {
                    for (row = 0; row < rows; row++)
                    {
                        for (col = 0; col < cols; col++)
                        {
                            for (layer = 0; layer < layers; layer++)
                            {
                                sw.Write(dtm[row, col]
                                    + "_" + layerthickness_m[row, col, layer]
                                    + "_" + bulkdensity[row, col, layer]
                                    + "_" + texture_kg[row, col, layer, 0]
                                    + "_" + texture_kg[row, col, layer, 1]
                                    + "_" + texture_kg[row, col, layer, 2]
                                    + "_" + texture_kg[row, col, layer, 3]
                                    + "_" + texture_kg[row, col, layer, 4]
                                    + "_" + young_SOM_kg[row, col, layer]
                                    + "_" + old_SOM_kg[row, col, layer]
                                    + ",");
                            }
                            sw.Write("\n");
                        }
                    }
                }
                catch { Debug.WriteLine(" issue with writing the content of the full output file for this timestep"); }
                sw.Close();
            }
        }

        public static double calc_slope_stdesc(int row_s, int col_s)
        {
            double slope_desc = 0, slope_temp = 0;
            if (dtm[row_s, col_s] != -9999)
            {
                for (int i = (-1); i <= 1; i++)
                {
                    for (int j = (-1); j <= 1; j++)
                    {
                        if (((row_s + i) >= 0) && ((row_s + i) < nr) && ((col_s + j) >= 0) && ((col_s + j) < nc) && !((i == 0) && (j == 0)))  //to stay within the grid and avoid the row col cell itself
                        {
                            if (dtm[row_s + i, col_s + j] != -9999) // if neighbour exists
                            {
                                if ((row_s != row_s + i) && (col_s != col_s + j)) { d_x = dx * Math.Sqrt(2); } else { d_x = dx; }
                                slope_temp = (dtm[row_s, col_s] - dtm[row_s + i, col_s + j]) / d_x;
                                if (slope_desc < slope_temp) { slope_desc = slope_temp; }
                            }
                        }
                    }
                }
            }

            slope_desc = Math.Atan(slope_desc); // slope in radians
            return (slope_desc);
        }

        public static void update_slope_and_aspect()
        {
            double slopemax, slope, slopetot;
            for (row = 0; row < nr; row++)
            {
                for (col = 0; col < nc; col++)
                {
                    if (dtm[row, col] != -9999)
                    {
                        slopemax = 0;
                        slope = 0;
                        slopetot = 0;

                        // Do slope analysis and Aspect Calculation first
                        if ((row - 1) >= 0)
                        {
                            if (dtm[row, col] > dtm[row - 1, col] && dtm[row - 1, col] != -9999) // North 0
                            {
                                slope = (dtm[row, col] - dtm[row - 1, col]) / dx;
                                if (slope > slopemax)
                                {
                                    slopemax = slope;
                                    slopetot++;
                                    aspect[row, col] = 0 * (3.141592654 / 180);
                                }
                            }
                        }

                        if ((row - 1) >= 0 & (col + 1) < nc)
                        {
                            if (dtm[row, col] > dtm[row - 1, col + 1] && dtm[row - 1, col + 1] != -9999) // Northeast 45
                            {
                                slope = (dtm[row, col] - dtm[row - 1, col + 1]) / (dx * Math.Sqrt(2));
                                if (slope > slopemax)
                                {
                                    slopemax = slope;
                                    slopetot++;
                                    aspect[row, col] = 45 * (3.141592654 / 180);
                                }
                            }
                        }

                        if ((col + 1) < nc)
                        {
                            if (dtm[row, col] > dtm[row, col + 1] && dtm[row, col + 1] != -9999) // East 90
                            {
                                slope = (dtm[row, col] - dtm[row, col + 1]) / dx;
                                if (slope > slopemax)
                                {
                                    slopemax = slope;
                                    slopetot++;
                                    aspect[row, col] = 90 * (3.141592654 / 180);
                                }
                            }
                        }
                        if ((row + 1) < nr & (col + 1) < nc)
                        {
                            if (dtm[row, col] > dtm[row + 1, col + 1] && dtm[row + 1, col + 1] != -9999) // SouthEast 135
                            {
                                slope = (dtm[row, col] - dtm[row + 1, col + 1]) / (dx * Math.Sqrt(2));
                                if (slope > slopemax)
                                {
                                    slopemax = slope;
                                    slopetot++;
                                    aspect[row, col] = 135 * (3.141592654 / 180);
                                }

                            }
                        }

                        if ((row + 1) < nr)
                        {
                            if (dtm[row, col] > dtm[row + 1, col] && dtm[row + 1, col] != -9999) // South 180
                            {
                                slope = (dtm[row, col] - dtm[row + 1, col]) / dx;
                                if (slope > slopemax)
                                {
                                    slopemax = slope;
                                    slopetot++;
                                    aspect[row, col] = 180 * (3.141592654 / 180);
                                }
                            }
                        }
                        if ((row + 1) < nr & (col - 1) >= 0)
                        {
                            if (dtm[row, col] > dtm[row + 1, col - 1] && dtm[row + 1, col - 1] != -9999) // SouthWest 225
                            {
                                slope = (dtm[row, col] - dtm[row + 1, col - 1]) / (dx * Math.Sqrt(2));
                                if (slope > slopemax)
                                {
                                    slopemax = slope;
                                    slopetot++;
                                    aspect[row, col] = 225 * (3.141592654 / 180);
                                }
                            }
                        }

                        if ((col - 1) >= 0)
                        {
                            if (dtm[row, col] > dtm[row, col - 1] && dtm[row, col - 1] != -9999) // West 270
                            {
                                slope = (dtm[row, col] - dtm[row, col - 1]) / dx;
                                if (slope > slopemax)
                                {
                                    slopemax = slope;
                                    slopetot++;
                                    aspect[row, col] = 270;
                                }
                            }
                        }

                        if ((row - 1) >= 0 & (col - 1) >= 0)
                        {
                            if (dtm[row, col] > dtm[row - 1, col - 1] && dtm[row - 1, col - 1] != -9999) // Northwest 315
                            {
                                slope = (dtm[row, col] - dtm[row - 1, col - 1]) / (dx * Math.Sqrt(2));
                                if (slope > slopemax)
                                {
                                    slopemax = slope;
                                    slopetot++;
                                    aspect[row, col] = 315 * (3.141592654 / 180);
                                }
                            }
                        }


                        if (slope > 0) slopeAnalysis[row, col] = slopemax;// Tom's: (slope/slopetot); ?
                        else { slopeAnalysis[row, col] = 0; }

                        // Convert slope to radians
                        slopeAnalysis[row, col] = System.Math.Atan(slopeAnalysis[row, col]);

                        //// test
                        //slopeAnalysis[row, col] = 0 * Math.PI / 180;
                        //aspect[row, col] = 0 * Math.PI / 180;

                    }
                }
            }
        }

    }
}
