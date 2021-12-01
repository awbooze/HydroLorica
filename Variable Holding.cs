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

        //if none of these fields get additional instructions, remove and rename the delegates
        //Ex: updateAllFields becomes UpdateAllFields
        protected void UpdateAllFields()
        {
            updateAllFields();
        }
        protected void UpdateStatusPannel()
        {
            updateStatusPannel();
        }
        protected void UpdateTimePannel()
        {
            updateTimePannel();
        }
        protected void UpdateFields()
        {
            updateVariables();
        }



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
    }
}
