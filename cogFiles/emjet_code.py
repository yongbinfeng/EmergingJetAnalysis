"""Generates code for EmergingJetAnalyzer"""
import cog
nonVectorDefault = '-1' # Default value for non vector variables
type_int         = 'int                     '
type_vint        = 'vector<int>             '
type_vvint       = 'vector< vector<int> >   '
type_float       = 'float                   '
type_vfloat      = 'vector<float>           '
type_vvfloat     = 'vector< vector<float> > '
event_vars = [
    ( "run                   " , type_int     ) ,
    ( "lumi                  " , type_int     ) ,
    ( "event                 " , type_int     ) ,
    ( "bx                    " , type_int     ) ,
    ( "met_pt                " , type_float   ) ,
    ( "met_phi               " , type_float   ) ,
]
jet_vars = [
    ( "jets_pt               " , type_vfloat  ) ,
    ( "jets_eta              " , type_vfloat  ) ,
    ( "jets_phi              " , type_vfloat  ) ,
    ( "jets_cef              " , type_vfloat  ) ,
    ( "jets_nef              " , type_vfloat  ) ,
    ( "jets_chf              " , type_vfloat  ) ,
    ( "jets_nhf              " , type_vfloat  ) ,
    ( "jets_phf              " , type_vfloat  ) ,
    ( "jets_nPromptTracks    " , type_vint    ) ,
    ( "jets_nDispTracks      " , type_vint    ) ,
    ( "jets_nSV              " , type_vint    ) ,
    ( "jets_medianLogIpSig   " , type_vfloat  ) ,
    ( "jets_missHits         " , type_vint    ) ,
    ( "jets_muonHits         " , type_vint    ) ,
    ( "jets_alphaMax         " , type_vfloat  ) ,
    ( "jets_nDarkPions       " , type_vint    ) ,
    ( "jets_minDRDarkPion    " , type_vfloat  ) ,
]
jet_track_vars = [
    ( "tracks_source         " , type_vvint   ) ,
    ( "tracks_pt             " , type_vvfloat ) ,
    ( "tracks_eta            " , type_vvfloat ) ,
    ( "tracks_phi            " , type_vvfloat ) ,
    ( "tracks_algo           " , type_vvint   ) ,
    ( "tracks_originalAlgo   " , type_vvint   ) ,
    ( "tracks_nHits          " , type_vvint   ) ,
    ( "tracks_nMissInnerHits " , type_vvint   ) ,
    ( "tracks_ipXY           " , type_vvfloat ) ,
    ( "tracks_ipZ            " , type_vvfloat ) , # Empty for now
    ( "tracks_ipXYSig        " , type_vvfloat ) ,
]
jet_vertex_vars = [
    ( "jet_vertex_source     " , type_vvint    ) ,
    ( "jet_vertex_x          " , type_vvfloat  ) ,
    ( "jet_vertex_y          " , type_vvfloat  ) ,
    ( "jet_vertex_z          " , type_vvfloat  ) ,
    ( "jet_vertex_xError     " , type_vvfloat  ) ,
    ( "jet_vertex_yError     " , type_vvfloat  ) ,
    ( "jet_vertex_zError     " , type_vvfloat  ) ,
    ( "jet_vertex_deltaR     " , type_vvfloat  ) ,
    ( "jet_vertex_Lxy        " , type_vvfloat  ) ,
    ( "jet_vertex_mass       " , type_vvfloat  ) ,
    ( "jet_vertex_chi2       " , type_vvfloat  ) ,
    ( "jet_vertex_ndof       " , type_vvfloat  ) ,
    ( "jet_vertex_pt2sum     " , type_vvfloat  ) ,
]
vertex_vars = [
    ( "vertex_source         " , type_vint    ) ,
    ( "vertex_x              " , type_vfloat  ) ,
    ( "vertex_y              " , type_vfloat  ) ,
    ( "vertex_z              " , type_vfloat  ) ,
    ( "vertex_xError         " , type_vfloat  ) ,
    ( "vertex_yError         " , type_vfloat  ) ,
    ( "vertex_zError         " , type_vfloat  ) ,
    ( "vertex_Lxy            " , type_vfloat  ) ,
    ( "vertex_mass           " , type_vfloat  ) ,
    ( "vertex_chi2           " , type_vfloat  ) ,
    ( "vertex_ndof           " , type_vfloat  ) ,
    ( "vertex_pt2sum         " , type_vfloat  ) ,
]
all_vars = event_vars + jet_vars + jet_track_vars + jet_vertex_vars + vertex_vars

def printLine(*args):
    line = "".join(args)
    # print line
    cog.outl(line)

def genOutputTree():
    # Output <typename> <varname>;
    for var in all_vars:
        varname = var[0]
        typename = var[1]
        printLine("%s %s;" % (typename, varname))

def genInit():
    # Output <varname>.clear();
    for var in all_vars:
        varname = var[0]
        typename = var[1]
        # Clear vectors
        if 'vector' in typename:
            printLine("%s.clear();" % varname)
        else:
            printLine("%s= %s;" % (varname, nonVectorDefault))

def genBranch():
    # printLine("#define BRANCH(tree, branch) (tree)->Branch(#branch, &branch);")
    # Output BRANCH(tree, <varname>);
    for var in all_vars:
        varname = var[0]
        typename = var[1]
        printLine("BRANCH(tree, %s);" % varname)
