class PtmConfig(object):
    def __init__(self):
        self.groups=[]
        self.release_distributions=[]

    def ptm_inp_text(self):
        sections=[self.global_information_text,
                  self.lines_text,
                  self.regions_text,
                  self.release_distributions_text,
                  self.release_timings_text,
                  self.behaviors_text,
                  self.output_text,
                  self.groups_text]
        blocks=["\n".join(lines) for lines in sections]
        return "\n".join(blocks)

    global_information_text=[
        "  GLOBAL INFORMATION",
        "   END_TIME = '2014-04-30 00:00:00'",
        "   RESTART_DIR = 'none'",
        "   TIME_STEP_SECONDS = 300.",
        "   -- deactivation logicals ---",
        "   REACHED_OPEN_BOUNDARY = 'true'",
        "   REACHED_FLOW_BOUNDARY = 'true'",
        "   ENTRAINED_BY_VOLUME_SINK = 'false'",
        "   CROSSED_LINE = 'false'",
        "   DEPOSITED_ON_BED = 'false'",
        "   CONSOLIDATED_ON_BED = 'false'",
        " ",
        "   -- kill logicals ---",
        "   REACHED_OPEN_BOUNDARY = 'false'",
        "   REACHED_FLOW_BOUNDARY = 'false'",
        "   ENTRAINED_BY_VOLUME_SINK = 'false'",
        "   CROSSED_LINE = 'false'",
        "   DEPOSITED_ON_BED = 'false'",
        "   CONSOLIDATED_ON_BED = 'false'",
        " "]
    lines_text=[
        "   -- line information --- ",
        "   NLINES = 0",
        "",
        " TRANSECT INFORMATION -- applies to tidal surfing",
        "   NTRANSECTS = 0",
        "",
        ]
    regions_text=[
        " REGION INFORMATION",
        "   NREGIONS = 1",
        "     -- -- region  1 --- ",
        "     -- REGION = 'mwt'",
        "     -- REGION_POLYGON_FILE = 'mwt.pol'",
        "     REGION = 'full_domain'",
        "     REGION_POLYGON_FILE = '../polygons/full_domain.pol'",
        "",
        ]

    release_distributions_text=[
        " RELEASE DISTRIBUTION INFORMATION",
        "   NRELEASE_DISTRIBUTION_SETS = 3",
        "",
        "   -- release distribution set 1 ---",
        "   RELEASE_DISTRIBUTION_SET = 'SAC' ",
        "   MIN_BED_ELEVATION_METERS = -99.",
        "   MAX_BED_ELEVATION_METERS =  99. ",
        "   HORIZONTAL_DISTRIBUTION = 'side'",
        "   SIDE_IDENTIFICATION_METHOD = 'find'",
        "     -- these must be within < 1 m.",
        "     XSTART = 619941.36",
        "     YSTART = 4293572.81",
        "     XEND   = 620029.72",
        "     YEND   = 4293754.71",
        "   NPARTICLE_ASSIGNMENT = 'flow'",
        "     DENSITY_SPECIFICATION_METHOD = 'constant'",
        "     PARTICLE_DENSITY = 0.0001",
        "     PARTICLE_DISTRIBUTION = 'random'",
        "   ZMIN_NON_DIM = 0.00",
        "   ZMAX_NON_DIM = 1.00",
        "   VERT_SPACING = 'uniform'",
        "",
        "   -- release distribution set 2 -- ",
        "   RELEASE_DISTRIBUTION_SET = 'SRV' ",
        "   MIN_BED_ELEVATION_METERS = -99.",
        "   MAX_BED_ELEVATION_METERS =  99. ",
        "   HORIZONTAL_DISTRIBUTION = 'side'",
        "   SIDE_IDENTIFICATION_METHOD = 'find'",
        "     -- these must be within < 1 m.",
        "     XSTART = 615232.95",
        "     YSTART = 4224588.29",
        "     XEND   = 615838.43",
        "     YEND   = 4224134.90",
        "   NPARTICLE_ASSIGNMENT = 'flow'",
        "     DENSITY_SPECIFICATION_METHOD = 'constant'",
        "     PARTICLE_DENSITY = 0.0001",
        "     PARTICLE_DISTRIBUTION = 'random'",
        "   ZMIN_NON_DIM = 0.00",
        "   ZMAX_NON_DIM = 1.00",
        "   VERT_SPACING = 'uniform'",
        "",
        "   -- release distribution set 3 ---",
        "   RELEASE_DISTRIBUTION_SET = 'PERCELL' ",
        "   MIN_BED_ELEVATION_METERS = -99.",
        "   MAX_BED_ELEVATION_METERS =  99. ",
        "   HORIZONTAL_DISTRIBUTION = 'region'",
        "   DISTRIBUTION_IN_REGION = 'cell'",
        "     -- RH not sure what this does, but distributed was not implemented.",
        "     -- so how about independent, the only other option.",
        "     CELL_RELEASE_TIMING = 'independent'",
        "     PARTICLE_NUMBER_CALCULATION_BASIS = 'volume'",
        "     -- aim high ",
        "     VOLUME_PER_PARTICLE_CUBIC_METERS = 10000.",
        "   ZMIN_NON_DIM = 0.00",
        "   ZMAX_NON_DIM = 1.00",
        "   VERT_SPACING = 'uniform'",
        ""]
    release_timings_text=[
        " RELEASE TIMING INFORMATION",
        "   NRELEASE_TIMING_SETS = 2",
        "",
        "   -- release timing set 1 ---        ",
        "   RELEASE_TIMING_SET = 'continuous'",
        "   INITIAL_RELEASE_TIME = '2014/04/02 00:00:00'",
        "   RELEASE_TIMING = 'continuous'",
        "   INACTIVATION_TIME = 'none'",
        "",
        "   -- release timing set 2 ---        ",
        "   RELEASE_TIMING_SET = 'once'",
        "   INITIAL_RELEASE_TIME = '2014/04/02 00:00:00'",
        "   RELEASE_TIMING = 'single'",
        "   INACTIVATION_TIME = 'none'",
        "",
        ]
    behaviors_text=[
        " -- BEHAVIOR INFORMATION",
        "   NBEHAVIOR_PROFILES = 0",
        "",
        "   NBEHAVIORS = 0",
        "",
        ]
    output_text=[
        " OUTPUT INFORMATION ",
        "   NOUTPUT_SETS = 1",
        "   -- output set 1 ---",
        "   OUTPUT_SET = '15min_output'",
        "   FLAG_LOG_LOGICAL = 'true'",
        "   BINARY_OUTPUT_INTERVAL_HOURS = 0.25",
        "   ASCII_OUTPUT_INTERVAL_HOURS = 'none' ",
        "   HISTOGRAM_OUTPUT_INTERVAL_HOURS = 'none'",
        "   STATISTICS_OUTPUT_INTERVAL_HOURS = 'none'",
        "   CONCENTRATION_OUTPUT_INTERVAL_HOURS = 'none'",
        "   REGION_COUNT_OUTPUT_INTERVAL_HOURS = 'none'",
        "   REGION_COUNT_UPDATE_INTERVAL_HOURS = 'none'",
        "   STATE_OUTPUT_INTERVAL_HOURS = 'none'",
        "",]
    @property
    def n_groups(self):
        return len(self.groups)

    @property
    def groups_text(self):
        pieces = [
            " PARTICLE GROUP INFORMATION ",
            " NGROUPS = %d"%self.n_groups
        ]
        for group in self.groups:
            pieces.extend(group)

        return pieces

        # "   -- group 1 ---",
        # "   GROUP = 'ryi'",
        # "   RELEASE_DISTRIBUTION_SET = 'SAC'",
        # "   RELEASE_TIMING_SET = 'continuous'",
        # "   PARTICLE_TYPE = 'none'",
        # "   BEHAVIOR_SET = 'none'",
        # "   OUTPUT_SET = '15min_output'",
        # "   OUTPUT_FILE_BASE = 'SAC'",
        # "",
        # "   -- group 2 ---",
        # "   GROUP = 'srv'",
        # "   RELEASE_DISTRIBUTION_SET = 'SRV'",
        # "   RELEASE_TIMING_SET = 'continuous'",
        # "   PARTICLE_TYPE = 'none'",
        # "   BEHAVIOR_SET = 'none'",
        # "   OUTPUT_SET = '15min_output'",
        # "   OUTPUT_FILE_BASE = 'SRV'",
        # "",
        # "   -- group 3 ---",
        # "   GROUP = 'initial'",
        # "   RELEASE_DISTRIBUTION_SET = 'PERCELL'",
        # "   REGION = 'full_domain'",
        # "   RELEASE_TIMING_SET = 'once'",
        # "   PARTICLE_TYPE = 'none'",
        # "   BEHAVIOR_SET = 'none'",
        # "   OUTPUT_SET = '15min_output'",
        # "   OUTPUT_FILE_BASE = 'INIT'",

##

from stompy.io.local import usgs_nwis
import numpy as np
from stompy.spatial import proj_utils
ll2utm=proj_utils.mapper('WGS84','EPSG:26910')

cache_dir="../particle_based/cache"

run_start=np.datetime64('2014-04-02')
run_stop =np.datetime64('2014-04-10')

stations=[
    dict(name='cli',code=11455315,desc="Cache Slough at S. Liberty Island")
]

for station in stations:
    ds=usgs_nwis.nwis_dataset(station['code'],run_start,run_stop,products=[10],
                              cache_dir=cache_dir,cache_only=True)
    if ds is None:
        print("No data for %s -- skip."%station['code'])
        continue
    meta=usgs_nwis.station_metadata(station['code'],cache_dir=cache_dir)
    ds['latitude']=(),meta['lat']
    ds['longitude']=(),meta['lon']
    ll=np.array( [ds.longitude.values,ds.latitude.values])
    ds['ll']=('xy',),ll
    xy=ll2utm(ll)
    ds['xy']=('xy',),xy
    ds.attrs['name']=station['name']
    ds.attrs['description']=station['desc']
    station['ds']=ds

##

pc=PtmConfig()

for station in stations:
    group=[
        "   GROUP = '%s'"%station['name'],
         "   RELEASE_DISTRIBUTION_SET = '%s'"%station['name'],
         "   RELEASE_TIMING_SET = 'continuous'",
         "   PARTICLE_TYPE = 'none'",
         "   BEHAVIOR_SET = 'none'",
         "   OUTPUT_SET = '15min_output'",
         "   OUTPUT_FILE_BASE = '%s'"%station['name'],
        ]
    pc.add_group(group)
