
from __future__ import print_function, division

import os, sys
import glob
import numpy as np
import datetime


from sqlalchemy import create_engine, Column, Integer, String, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

## import from local files
## Boilerplate path hack to give access to full clustered_SNe package
import sys
if __package__ is None:
    if os.pardir not in sys.path[0]:
        file_dir = os.getcwd()
        sys.path.insert(0, os.path.join(file_dir, 
                                        os.pardir, 
                                        os.pardir))

from clustered_SNe.analysis.constants import m_proton, pc, yr, M_solar, \
                                    metallicity_solar
    
from clustered_SNe.analysis.parse import parse_run, Overview, Inputs, \
                                        extract_masses_momenta_raw


#####################
# Required objects, part 1
# hardcodes database name -- probably okay?
# check later: is this the best way to do this?

import warnings
warnings.warn("`session' from database_helpers can only write using 1 process at a time",
    UserWarning)

engine = create_engine('sqlite:///clustered_SNe.db', echo=False)
Base = declarative_base()


######### TABLES ############
class Simulation(Base):
    __tablename__ = "simulations"
    id          = Column(String, primary_key = True)
    data_dir    = Column(String)
    
    metallicity             = Column(Float)
    background_density      = Column(Float)
    background_temperature  = Column(Float)
    with_cooling            = Column(Integer)
    cooling_type            = Column(String)
    num_SNe                 = Column(Integer)
    cluster_mass            = Column(Float)
    seed                    = Column(Integer)
    mass_loss               = Column(String)
    
    E_R_kin                 = Column(Float)
    E_R_int                 = Column(Float)
    M_R                     = Column(Float)
    R                       = Column(Float)
    t                       = Column(Float)
    momentum                = Column(Float)
    num_checkpoints         = Column(Integer)
    
    last_updated            = Column(String)
    
    def __repr__(self):
        return "<Simulation: {0}, last updated: {1}".format(self.id,
                                                            self.last_updated)

    @classmethod
    def from_last_run(cls, data_dir, last_run):
        id = last_run.id

        metallicity             = last_run.overview.metallicity
        background_density      = last_run.overview.background_density
        background_temperature  = last_run.overview.background_temperature
        with_cooling            = last_run.overview.with_cooling
        cooling_type            = last_run.overview.cooling_type
        num_SNe                 = last_run.overview.num_SNe
        cluster_mass            = last_run.overview.cluster_mass
        seed                    = last_run.overview.seed
        mass_loss               = last_run.overview.mass_loss
        
        if num_SNe > 1:
            extraction_index = np.argmax(last_run.momentum)

            E_R_kin                 = last_run.E_R_kin[ extraction_index]
            E_R_int                 = last_run.E_R_int[ extraction_index]
            M_R                     = last_run.M_R[     extraction_index]
            R                       = last_run.R_shock[ extraction_index]

            t                       = last_run.times[   extraction_index]
            momentum                = last_run.momentum[extraction_index]

        else:
            E_R_kin                 = 0 
            E_R_int                 = 0 
            M_R                     = 0 
            R                       = 0 

            t                       = 0 
            momentum                = 0             

        num_checkpoints = len(last_run.filenames)

        last_updated    = str(datetime.datetime.now())

        simulation = cls(id=id, data_dir=data_dir,
            metallicity             = metallicity,
            background_density      = background_density,
            background_temperature  = background_temperature,
            with_cooling            = with_cooling,
            cooling_type            = cooling_type,
            num_SNe                 = num_SNe,
            cluster_mass            = cluster_mass, 
            seed                    = seed,
            mass_loss               = mass_loss, 
            E_R_kin                 = E_R_kin,
            E_R_int                 = E_R_int,
            M_R                     = M_R,
            R                       = R,
            t                       = t,
            momentum                = momentum,
            num_checkpoints         = num_checkpoints,
            last_updated            = last_updated)
        return simulation

    def add_to_table(self):
        session.add(self)

    def update_to_table(self):
        # Only updates the things I *think* are going to change
        # things like cluster mass *shouldn't* change, 
        # but if you did change that, things could break
        session.query(self.__class__).\
            filter(self.__class__.id==self.id).\
            update({
                "E_R_kin": self.E_R_kin,
                "E_R_int": self.E_R_int,
                "M_R":     self.M_R,
                "R":       self.R,
                "t":       self.t,
                "momentum": self.momentum,
                "num_checkpoints": self.num_checkpoints,
                "last_updated": self.last_updated
            })

    def add_or_update_to_table(self):
        existing_entry = session.query(self.__class__).\
            get(self.id)

        if existing_entry is None:
            self.add_to_table()
        else:
            self.update_to_table()



class Simulation_Inputs(Base):
    __tablename__ = "inputs"
    id = Column(String, primary_key = True)
    
    T_Start             = Column(Float)
    T_End               = Column(Float)
    Num_Reports         = Column(Integer)
    Num_Checkpoints     = Column(Integer)
    Use_Logtime         = Column(Integer)
    
    Num_R               = Column(Integer)
    R_Min               = Column(Float)
    R_Max               = Column(Float)
    Log_Zoning          = Column(Integer)
    Log_Radius          = Column(Integer)
    
    CFL                 = Column(Integer)
    PLM                 = Column(Integer)
    RK2                 = Column(Integer)
    H_0                 = Column(Float)
    H_1                 = Column(Float)
    Riemann_Solver      = Column(Integer)
    Density_Floor       = Column(Float)
    Pressure_Floor      = Column(Float)
    
    With_Cooling        = Column(Integer)
    Cooling_Type        = Column(String)
    Adiabatic_Index     = Column(Float)
    
    ICs                 = Column(String)
    mass_loss           = Column(String)
    
    def __repr__(self):
        return "<{0} for id: {1}>".format(self.__class__.__name__,
                                          self.id)
    
    @staticmethod
    def from_Inputs(id, inputs):
        
        simulation_inputs = Simulation_Inputs(id=id,
            T_Start         = inputs.T_Start,
            T_End           = inputs.T_End,
            Num_Reports     = inputs.Num_Reports,
            Num_Checkpoints = inputs.Num_Checkpoints,
            Use_Logtime     = inputs.Use_Logtime,

            Num_R           = inputs.Num_R,
            R_Min           = inputs.R_Min,
            R_Max           = inputs.R_Max,
            Log_Zoning      = inputs.Log_Zoning,
            Log_Radius      = inputs.Log_Radius,

            CFL             = inputs.CFL,
            PLM             = inputs.PLM,
            RK2             = inputs.RK2,
            H_0             = inputs.H_0,
            H_1             = inputs.H_1,
            Riemann_Solver  = inputs.Riemann_Solver,
            Density_Floor   = inputs.Density_Floor,
            Pressure_Floor  = inputs.Pressure_Floor,

            With_Cooling    = inputs.With_Cooling,
            Cooling_Type    = inputs.Cooling_Type,

            Adiabatic_Index = inputs.Adiabatic_Index,

            ICs             = inputs.ICs,

            mass_loss       = inputs.mass_loss,
        )
        return simulation_inputs

    def add_to_table(self):
        session.add(self)

    def update_to_table(self):
        # none of the inputs should change;
        # if they do, you're probably breaking other things
        pass

    def add_or_update_to_table(self):
        existing_entry = session.query(self.__class__).\
            get(self.id)

        if existing_entry is None:
            self.add_to_table()
        else:
            self.update_to_table()


class Simulation_Status(Base):
    __tablename__ = "status"

    id = Column(String, primary_key=True)
    data_dir = Column(String)
    status   = Column(String)

    # It'd be really nice to enforce these options
    possible_statuses = set(["Complete",
                             "Ready",
                             "Error",
                             "Running",
                             "Unknown"])
    # Complete = finished AND converged
    # Ready = Can be restarted, but isn't restarted yet
    #         - might not NEED to be restarted, but could be
    # Error = failed -- needs to be reinitialized and deleted
    # Running = restarted (or in queue...)
    # Unknown = everything else...


    # do all these static/class methods probably need to be static/class methods
    @staticmethod
    def is_converged(id, data_dir):
        overview = Overview(os.path.join(data_dir, id.rstrip("_") + "_overview.dat"))
        if overview.num_SNe == 0:
            return True

        last_run = parse_run(data_dir, id)

        momenta = last_run.momentum
        if (momenta.size < 25):
            earlier_momentum = momenta[0]
        else:
            earlier_momentum = momenta[-25]
        current_momentum = momenta[-1]
        tolerance = .01
        if (current_momentum > ((1+tolerance) * earlier_momentum)):
            return False
        else:
            return True

    @classmethod
    def mappable_is_converged(cls, id_and_data_dir):
        id, data_dir = id_and_data_dir
        return cls.is_converged(id, data_dir)


    @staticmethod
    def is_last_checkpoint_x99(id, data_dir):
        checkpoints = glob.glob(os.path.join(data_dir, id + "_checkpoint_*.dat"))
        checkpoints = sorted(checkpoints)

        last_checkpoint = checkpoints[-1]
        last_checkpoint_num = int(last_checkpoint.split("_")[-1].strip(".dat")) % 100
        if last_checkpoint_num != 99:
            return False
        else:
            return True

    @classmethod
    def mappable_is_last_checkpoint_x99(cls, id_and_data_dir):
        id, data_dir = id_and_data_dir
        return cls.is_last_checkpoint_x99(id, data_dir)


    def create_new_batch_file(self, 
                              scripts_dir     = "../scripts/",
                              destination_dir = "../scripts/",
                              verbose=False):
        with open(os.path.join(scripts_dir, "restart.batch")) as f:
            text = f.read()
        with open(os.path.join(destination_dir, "restart.batch."+self.id), mode="w") as f2:
            text2 = text[:-1] + self.id
            if verbose:
                print(text2)
            f2.write(text2)

    def add_to_table(self):
        session.add(self)

    def update_to_table(self):
        session.query(self.__class__).\
            filter(self.__class__.id==self.id).\
            update({"status":self.status})

    def add_or_update_to_table(self):
        assert(self.status in self.possible_statuses)

        existing_entry = session.query(self.__class__).\
            get(self.id)

        if existing_entry is None:
            self.add_to_table()
        elif existing_entry.status == "Complete":
            return
        elif existing_entry.status == "Error":
            return
        elif (existing_entry.status == "Running") and (self.status == "Unknown"):
            return
        else:
            self.update_to_table()

    def switch_to_running(self):
        self.status = "Running"
        self.add_or_update_to_table()

    def create_new_batch_file_and_switch_to_running(self,
        scripts_dir     = "../scripts/",
        destination_dir = "../scripts/",
        verbose = False):

        self.create_new_batch_file(scripts_dir = scripts_dir,
            destination_dir = destination_dir,
            verbose = verbose)

        self.switch_to_running()

    def __repr__(self):
        return "<{0} for id: {1}>".format(self.__class__.__name__,
                                          self.id)




#####################
# Required objects, part 2
# check later: is this the best way to do this?
Base.metadata.create_all(engine)

Session = sessionmaker(bind=engine)
session = Session()



####################

def extract_masses_momenta(density, metallicity):
    """Fetches the cluster mass and max momentum from SQL table.
    For more options, check out extract_masses_momenta_raw from `parse.py`"""

    ids = []
    masses = []
    momenta = []
    print("density: ", density)
    print("metallicity: ", metallicity)

    for simulation in session.query(Simulation).\
                    filter(Simulation.metallicity==metallicity):
        if np.isclose(simulation.background_density, density, rtol=.1, atol=0):
            ids.append(simulation.id)
            masses.append(simulation.cluster_mass)
            momenta.append(simulation.momentum)
    
    if len(ids) == 0:
        print("No matching files found!")
        return
    
    ids     = np.array(ids)
    masses  = np.array(masses)
    momenta = np.array(momenta)


    sorted_indices = np.argsort(masses)
    masses  = masses[  sorted_indices]
    momenta = momenta[ sorted_indices]
    ids     = ids[     sorted_indices]
    return masses, momenta, ids




