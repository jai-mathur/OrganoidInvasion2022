from cc3d.core.PySteppables import *
import numpy as np

front_angle = 3.14/4
F_ECM = 1.0
thresh = 0.05 # 5% threshold for the hill function
f_R = {{f_R}}
P0 = {{P0}}#This one is used with the active ECM feedback mechanism
PL = {{PL}}#This is directly feeding it values.
PF = 1.0
class CalculateARandFrontLeaders(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)
        
        
    def start(self):
        """
        any code in the start function runs before MCS=0
        """
        

    def step(self,mcs):
        """
        type here the code that will run every frequency MCS
        :param mcs: current Monte Carlo step
        """
        x_min = 0.0
        x_max = 0.0
        y_min = 0.0
        y_max = 0.0
        z_min = 0.0
        z_max = 0.0
        
        org_com = np.array([0.0, 0.0, 0.0])
        fol_com = np.array([0.0, 0.0, 0.0])
        leader_com = np.array([0.0, 0.0, 0.0])
        
        front_count = 0
        
        for cell in self.cell_list:
            x_min = (x_min < cell.xCOM)*(x_min - cell.xCOM) + cell.xCOM
            y_min = (y_min < cell.yCOM)*(y_min - cell.yCOM) + cell.yCOM
            z_min = (z_min < cell.zCOM)*(z_min - cell.zCOM) + cell.zCOM
            
            z_max = (z_max > cell.zCOM)*(z_max - cell.zCOM) + cell.zCOM
            y_max = (y_max > cell.yCOM)*(y_max - cell.yCOM) + cell.yCOM
            x_max = (x_max > cell.xCOM)*(x_max - cell.xCOM) + cell.xCOM
            
            org_com += np.array([cell.xCOM, cell.yCOM, cell.zCOM])
            
            if cell.type == self.cell_type.FOLLOWER:
                fol_com += np.array([cell.xCOM, cell.yCOM, cell.zCOM])
                
            if cell.type == self.cell_type.LEADER:
                leader_com += np.array([cell.xCOM, cell.yCOM, cell.zCOM])            
                
        org_com = org_com/len(self.cell_list)
        fol_com = fol_com/len(self.cell_list_by_type(self.cell_type.FOLLOWER))
        leader_com = leader_com/len(self.cell_list_by_type(self.cell_type.LEADER))
        
        self.shared_steppable_vars['current_org_com'] = org_com
        self.shared_steppable_vars['current_fol_com'] = fol_com
        self.shared_steppable_vars['current_leader_com'] = leader_com
        
        self.shared_steppable_vars['org_AR'] = ( ( (y_max-y_min) + (z_max-z_min) )/2 )/(x_max-x_min) # the num and denom both had a divide by 2 that got cancelled out

        for cell in self.cell_list_by_type(self.cell_type.LEADER):
            vec = -fol_com + np.array([cell.xCOM, cell.yCOM, cell.zCOM])
            x_unit_vector = np.array([1,0,0])
            dot_prod = np.dot(vec,x_unit_vector)
            cos_theta = dot_prod/self.vector_norm(vec)
            
            if (dot_prod > 40*np.cos(front_angle) ):#40*cos(theta) defines our conical region
                front_count += 1
          
        self.shared_steppable_vars['front_leaders'] = front_count



class Calculate_P(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)
        
    def start(self):
        """
        any code in the start function runs before MCS=0
        """
        initial_org_com = np.array([0.0, 0.0, 0.0])
        initial_L_com = np.array([0.0, 0.0, 0.0])
        initial_fol_com = np.array([0.0, 0.0, 0.0])
        self.shared_steppable_vars['initial_org_com'] = np.array([0.0, 0.0, 0.0])
        self.shared_steppable_vars['initial_leader_com'] = np.array([0.0, 0.0, 0.0])
        self.shared_steppable_vars['initial_fol_com'] = np.array([0.0, 0.0, 0.0])

        for cell in self.cell_list:
            initial_org_com += np.array([cell.xCOM, cell.yCOM, cell.zCOM])
        initial_org_com /= len(self.cell_list)
        self.shared_steppable_vars['initial_org_com'] = initial_org_com
        
        for cell in self.cell_list_by_type(self.cell_type.LEADER):
            vec = np.array([cell.xCOM,cell.yCOM,cell.zCOM]) - initial_org_com
            if self.vector_norm(vec) > 23:
                cell.type = self.cell_type.FOLLOWER

        for cell in self.cell_list_by_type(self.cell_type.LEADER):
            initial_L_com += np.array([cell.xCOM, cell.yCOM, cell.zCOM])
        initial_L_com /= len(self.cell_list_by_type(self.cell_type.LEADER))
        self.shared_steppable_vars['initial_leader_com'] = initial_L_com

        for cell in self.cell_list_by_type(self.cell_type.FOLLOWER):
            initial_fol_com += np.array([cell.xCOM, cell.yCOM, cell.zCOM])
        initial_fol_com /= len(self.cell_list_by_type(self.cell_type.FOLLOWER))
        self.shared_steppable_vars['initial_fol_com'] = initial_fol_com
                        

        print ("percentage of leaders", 100*len(self.cell_list_by_type(self.cell_type.LEADER) )/len(self.cell_list) )
        print ("num of leaders", len(self.cell_list_by_type(self.cell_type.LEADER) )  )

        for cell in self.cell_list:

            
            if cell.type == self.cell_type.LEADER:
                # cell.dict['P0'] = 10.0
                cell.dict['beta'] = 0.0
                
            if cell.type == self.cell_type.FOLLOWER:
                cell.dict['P'] = 0.0
##########################################################################################################################

        self.plot_win1 = self.add_new_plot_window(
            title='Average org speed',
            x_axis_title='MonteCarlo Step (MCS)',
            y_axis_title='org speed',
            x_scale_type='linear',
            y_scale_type='linear',
            grid=True # only in 3.7.6 or higher
        )
        self.plot_win1.add_plot("OrgSpeed", style='Dots', color='red', size=5)                

        self.plot_win2 = self.add_new_plot_window(
            title='Average Leader Speed',
            x_axis_title='MonteCarlo Step (MCS)',
            y_axis_title='Leader speed',
            x_scale_type='linear',
            y_scale_type='linear',
            grid=True # only in 3.7.6 or higher
        )
        self.plot_win2.add_plot("LeaderSpeed", style='Dots', color='red', size=5)                

        self.plot_win3 = self.add_new_plot_window(
            title='Average Follower Speed',
            x_axis_title='MonteCarlo Step (MCS)',
            y_axis_title='Follower Speed',
            x_scale_type='linear',
            y_scale_type='linear',
            grid=True # only in 3.7.6 or higher
        )
        self.plot_win3.add_plot("FolSpeed", style='Dots', color='red', size=5)                

        self.plot_win4 = self.add_new_plot_window(
            title='Organoid AR',
            x_axis_title='MonteCarlo Step (MCS)',
            y_axis_title='Org AR',
            x_scale_type='linear',
            y_scale_type='linear',
            grid=True # only in 3.7.6 or higher
        )
        self.plot_win4.add_plot("OrgAR", style='Dots', color='red', size=5)                


        self.plot_win5 = self.add_new_plot_window(
            title='Leader cell localization',
            x_axis_title='MonteCarlo Step (MCS)',
            y_axis_title='% of leader cell in front',
            x_scale_type='linear',
            y_scale_type='linear',
            grid=True # only in 3.7.6 or higher
        )
        self.plot_win5.add_plot("FrontLeadersPercent", style='Dots', color='red', size=5)                 

                    

    def step(self,mcs):
        """
        type here the code that will run every frequency MCS
        :param mcs: current Monte Carlo step
        """
        # CCCP = 1.0 Cell-Cell Communication Parameter
        # delta_t = 1/60*90/1 1 MCS = 1 minute
        global P0
        global f_R
        tot_leaders = len(self.cell_list_by_type(self.cell_type.LEADER) )

        for cell in self.cell_list_by_type(self.cell_type.LEADER):
            free_boundary_norm = 0.0
            flag_boundary = 0
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if not neighbor:
                    fbn = common_surface_area/cell.surface #fbn =  free boundary normalized
                    cell.dict['beta'] = fbn**4/(fbn**4 + thresh**4)
                    flag_boundary = 1
            if flag_boundary == 0:
                cell.dict['beta'] = 0.0

        n_activated_leader = 0

        for cell in self.cell_list_by_type(self.cell_type.LEADER):
            if cell.dict['beta'] > 0.95:
                n_activated_leader += 1
        fr_net = 0.0
        for cell in self.cell_list_by_type(self.cell_type.LEADER):
            beta = cell.dict['beta']
            fr_net = n_activated_leader/len(self.cell_list)*beta*f_R
            if beta == 0.0:
                P = P0 
            else:
                P = P0*(1.0 - 2.71**(-fr_net/F_ECM) )
            cell.dict['P'] = P 
            
            

        PF = 0.0
        beta_tot = 0.0
        for cell in self.cell_list_by_type(self.cell_type.LEADER):
            beta = cell.dict['beta']
            PF += cell.dict['P']*beta
            beta_tot += beta
        if beta_tot != 0.0:
            PF = PF/beta_tot

        for cell in self.cell_list_by_type(self.cell_type.FOLLOWER):
            cell.dict['P'] = PF


        # #print("266")
        if mcs == 0:        
            # self.shared_steppable_vars['current_LE_x'] = LE_x

            self.plot_win4.add_data_point('OrgAR', mcs, self.shared_steppable_vars['org_AR'])
            self.plot_win5.add_data_point('FrontLeadersPercent', mcs, self.shared_steppable_vars['front_leaders']/tot_leaders)

            #print("269")
        elif mcs%50 == 0:        
            delta_org = self.shared_steppable_vars['current_org_com'] - self.shared_steppable_vars['initial_org_com'] # NOT INSTANTANEOUS, BUT AVG
            delta_leader = self.shared_steppable_vars['current_leader_com'] - self.shared_steppable_vars['initial_leader_com'] # NOT INSTANTANEOUS, BUT AVG
            delta_fol = self.shared_steppable_vars['current_fol_com'] - self.shared_steppable_vars['initial_fol_com'] # NOT INSTANTANEOUS, BUT AVG
            
            self.plot_win1.add_data_point('OrgSpeed', mcs, (delta_org[0])/mcs)
            self.plot_win2.add_data_point('LeaderSpeed', mcs, (delta_leader[0])/mcs)
            self.plot_win3.add_data_point('FolSpeed', mcs, (delta_fol[0])/mcs)
            self.plot_win4.add_data_point('OrgAR',mcs, self.shared_steppable_vars['org_AR'])
            self.plot_win5.add_data_point('FrontLeadersPercent',mcs, self.shared_steppable_vars['front_leaders']/tot_leaders)


    def finish(self):
        """
        Finish Function is called after the last MCS
        """
        if self.output_dir is not None:
            png_output_path = Path(self.output_dir).joinpath("Plot_org_speed"+ ".png")
            self.plot_win1.save_plot_as_png(png_output_path)
            csv_output_path = Path(self.output_dir).joinpath("data_org_speed"+ ".csv")
            self.plot_win1.save_plot_as_data(csv_output_path, CSV_FORMAT)

            png_output_path = Path(self.output_dir).joinpath("Plot_leader_speed"+ ".png")
            self.plot_win2.save_plot_as_png(png_output_path)
            csv_output_path = Path(self.output_dir).joinpath("data_leader_speed"+ ".csv")
            self.plot_win2.save_plot_as_data(csv_output_path, CSV_FORMAT)

            png_output_path = Path(self.output_dir).joinpath("Plot_fol_speed"+ ".png")
            self.plot_win3.save_plot_as_png(png_output_path)
            csv_output_path = Path(self.output_dir).joinpath("data_fol_speed"+ ".csv")
            self.plot_win3.save_plot_as_data(csv_output_path, CSV_FORMAT)

            png_output_path = Path(self.output_dir).joinpath("Plot_orgAR"+ ".png")
            self.plot_win4.save_plot_as_png(png_output_path)
            csv_output_path = Path(self.output_dir).joinpath("data_orgAR"+ ".csv")
            self.plot_win4.save_plot_as_data(csv_output_path, CSV_FORMAT)

            png_output_path = Path(self.output_dir).joinpath("Plot_front_leaders"+ ".png")
            self.plot_win5.save_plot_as_png(png_output_path)
            csv_output_path = Path(self.output_dir).joinpath("data_front_leaders"+ ".csv")
            self.plot_win5.save_plot_as_data(csv_output_path, CSV_FORMAT)



class Spheroid_invasionSteppable(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)

    def start(self):
        """
        any code in the start function runs before MCS=0
        """
#         initial_pos_population = np.array([0.0,0.0,0.0])
        for cell in self.cell_list:
            cell.targetVolume = cell.volume
            cell.lambdaVolume = 1.0
            cell.targetSurface = 1.1*cell.surface
            cell.lambdaSurface = 0.1

    def step(self,mcs):
        """
        type here the code that will run every frequency MCS
        :param mcs: current Monte Carlo step
        """
###########################################################################################################################        
#        global PL# these are to be used for simulations without active ECM feedback
                 # WHEN YOU WANT FEEDBACK, ACTIVATE THE DICT['P'] STATEMENTS FOR P TO BE READ THROUGH THE MECHANISM 
                 # ELSE USE THE -1.0*PL & -1.0 FOR FOLLOWERS STATEMENT FOR DIRECT ASSIGNMENT OF  PL AND PF.
###########################################################################################################################        
        for cell in self.cell_list_by_type(self.cell_type.LEADER):

            cell.lambdaVecX = -1.0*cell.dict['P']
            # cell.lambdaVecX = -1.0*PL
            cell.lambdaVecY = 0.0
            cell.lambdaVecZ = 0.0
        
        for cell in self.cell_list_by_type(self.cell_type.FOLLOWER):

            cell.lambdaVecX = -1*( (cell.dict['P'] <=1)*(1 - cell.dict['P']) + cell.dict['P'] )
            # cell.lambdaVecX = -1.0
            cell.lambdaVecY = 0.0
            cell.lambdaVecZ = 0.0

        
        
        
        


    def finish(self):
        """
        Finish Function is called after the last MCS
        """