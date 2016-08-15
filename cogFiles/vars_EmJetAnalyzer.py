"""Generates code for EmergingJetAnalyzer"""

############################################################
# Cog file boiler plate
# Use outputline() instead of cog.outl()
standalone = False
if __name__=='__main__':
    print "Running standalone mode - Printing to screen"
    standalone = True
else:
    try: import cog
    except: ImportError

def outputline(line):
    if not standalone:
        cog.outl(line)
    else:
        print(line)
# End Cog file boiler plate
############################################################

# TYPESTRLENGTH = 30

def pad(input_string, width=20):
    '''Pads input string with spaces to given width'''
    format_string = '{0: <%d}' % width
    output_string = format_string.format(input_string)
    return output_string

non_vector_default = '-1' # Default value for non vector variables
from collections import namedtuple
# Var.level defines the vector level when storing in TTree
# e.g. cpptype="int", level=2 --> "vector<vector<int > >"
Var = namedtuple('Var', ['name', 'cpptype', 'level', ] )
event_vars = [
    Var("run"                 , "int"   , 0 , ) ,
    Var("lumi"                , "int"   , 0 , ) ,
    Var("event"               , "int"   , 0 , ) ,
    Var("bx"                  , "int"   , 0 , ) ,
    Var("nVtx"                , "int"   , 0 , ) ,
    Var("nGoodVtx"            , "int"   , 0 , ) ,
    Var("nTrueInt"            , "int"   , 0 , ) ,
    Var("met_pt"              , "float" , 0 , ) ,
    Var("met_phi"             , "float" , 0 , ) ,
]
jet_vars = [
    Var("index"               , "int"   , 1 , ) ,
    Var("source"              , "int"   , 1 , ) ,
    Var("pt"                  , "float" , 1 , ) ,
    Var("eta"                 , "float" , 1 , ) ,
    Var("phi"                 , "float" , 1 , ) ,
    Var("cef"                 , "float" , 1 , ) ,
    Var("nef"                 , "float" , 1 , ) ,
    Var("chf"                 , "float" , 1 , ) ,
    Var("nhf"                 , "float" , 1 , ) ,
    Var("phf"                 , "float" , 1 , ) ,
    # Var("nPromptTracks"       , "int"   , 1 , ) ,
    # Var("nDispTracks"         , "int"   , 1 , ) ,
    # Var("nSV"                 , "int"   , 1 , ) ,
    # Var("medianLogIpSig"      , "float" , 1 , ) ,
    Var("missHits"            , "int"   , 1 , ) ,
    Var("muonHits"            , "int"   , 1 , ) ,
    Var("alphaMax"            , "float" , 1 , ) ,
    Var("nDarkPions"          , "int"   , 1 , ) ,
    Var("nDarkGluons"          , "int"   , 1 , ) ,
    Var("minDRDarkPion"       , "float" , 1 , ) ,
]
jet_track_vars = [
    Var("index"               , "int"   , 2 , ) ,
    Var("source"              , "int"   , 2 , ) ,
    Var("jet_index"           , "int"   , 2 , ) ,
    Var("vertex_index"        , "int"   , 2 , ) ,
    Var("vertex_weight"       , "float" , 2 , ) ,
    Var("nHitsInFrontOfVert"  , "int"   , 2 , ) ,
    Var("missHitsAfterVert"   , "int"   , 2 , ) ,
    Var("pt"                  , "float" , 2 , ) ,
    Var("eta"                 , "float" , 2 , ) ,
    Var("phi"                 , "float" , 2 , ) ,
    Var("pca_r"               , "float" , 2 , ) ,
    Var("pca_eta"             , "float" , 2 , ) ,
    Var("pca_phi"             , "float" , 2 , ) ,
    Var("quality"             , "int"   , 2 , ) ,
    Var("algo"                , "int"   , 2 , ) ,
    Var("originalAlgo"        , "int"   , 2 , ) ,
    Var("nHits"               , "int"   , 2 , ) ,
    Var("nMissInnerHits"      , "int"   , 2 , ) ,
    Var("nTrkLayers"          , "int"   , 2 , ) ,
    Var("nMissInnerTrkLayers" , "int"   , 2 , ) ,
    Var("nMissOuterTrkLayers" , "int"   , 2 , ) ,
    Var("nMissTrkLayers"      , "int"   , 2 , ) ,
    Var("nPxlLayers"          , "int"   , 2 , ) ,
    Var("nMissInnerPxlLayers" , "int"   , 2 , ) ,
    Var("nMissOuterPxlLayers" , "int"   , 2 , ) ,
    Var("nMissPxlLayers"      , "int"   , 2 , ) ,
    Var("ipXY"                , "float" , 2 , ) ,
    Var("ipZ"                 , "float" , 2 , ) , # Empty for now
    Var("ipXYSig"             , "float" , 2 , ) ,
    Var("dRToJetAxis"         , "float" , 2 , ) ,
    Var("distanceToJet"       , "float" , 2 , ) ,
]
jet_vertex_vars = [
    Var("index"               , "int"   , 2 , ) ,
    Var("source"              , "int"   , 2 , ) ,
    Var("jet_index"           , "int"   , 2 , ) ,
    Var("x"                   , "float" , 2 , ) ,
    Var("y"                   , "float" , 2 , ) ,
    Var("z"                   , "float" , 2 , ) ,
    Var("xError"              , "float" , 2 , ) ,
    Var("yError"              , "float" , 2 , ) ,
    Var("zError"              , "float" , 2 , ) ,
    Var("deltaR"              , "float" , 2 , ) ,
    Var("Lxy"                 , "float" , 2 , ) ,
    Var("mass"                , "float" , 2 , ) ,
    Var("chi2"                , "float" , 2 , ) ,
    Var("ndof"                , "float" , 2 , ) ,
    Var("pt2sum"              , "float" , 2 , ) ,
]
all_vars = event_vars + jet_vars + jet_track_vars + jet_vertex_vars
# Pad Var fields with appropriate number of spaces
namelength_list = [len(var.name) for var in all_vars]
maxnamelength = max(namelength_list)

def var_to_dict(var):
    vardict = var._asdict()
    # Populate vardict fields, padding with appropriate number of spaces
    vardict['name']       = pad(var.name       , maxnamelength+1 )
    vardict['cpptype']    = pad(var.cpptype    , 5+1             )
    if   var.level == 0: branchtype = '%s'                  % var.cpptype
    elif var.level == 1: branchtype = 'vector<%s>'          % var.cpptype
    elif var.level == 2: branchtype = 'vector<vector<%s> >' % var.cpptype
    vardict['branchtype'] = pad(branchtype , 5+18          )
    return vardict

def make_fullname_builder(prefix="", postfix=""):
    """Return lambda that adds input['fullname'] = prefix + input['name'] to input dictionary"""
    # Hacky way to make a lambda that returns the dictionary after updating it
    return lambda x : x.update({'fullname': prefix + x['name'] + postfix}) or x

# Turn list of Var objects, to dictionary objects
event_vardicts      = map( var_to_dict, event_vars      )
jet_vardicts        = map( var_to_dict, jet_vars        )
jet_track_vardicts  = map( var_to_dict, jet_track_vars  )
jet_vertex_vardicts = map( var_to_dict, jet_vertex_vars )
# Add appropriate prefix for variable names
event_vardicts      = map( make_fullname_builder(""        ) , event_vardicts      )
jet_vardicts        = map( make_fullname_builder("jet_"    ) , jet_vardicts        )
jet_track_vardicts  = map( make_fullname_builder("track_"  ) , jet_track_vardicts  )
jet_vertex_vardicts = map( make_fullname_builder("vertex_" ) , jet_vertex_vardicts )
# for vardict in event_vardicts      : vardict['prefix'] = ""
# for vardict in jet_vardicts        : vardict['prefix'] = "jet_"
# for vardict in jet_track_vardicts  : vardict['prefix'] = "track_"
# for vardict in jet_vertex_vardicts : vardict['prefix'] = "vertex_"
all_vardicts = event_vardicts + jet_vardicts + jet_track_vardicts + jet_vertex_vardicts

from string import Template
def replaceSingleLine(template_string, vardict):
    t = Template(template_string)
    tsub = t.substitute(vardict)
    outputline(tsub)

def gen_OutputTree():
    """Generate OutputTree class declaration for OutputTree.h"""
    # Output <typename> <varname>;
    for vardict in all_vardicts:
        varname = vardict['fullname']
        typename = vardict['branchtype']
        outputline("%s %s;" % (typename, varname))

def gen_Init():
    """Generate Init() for OutputTree.h"""
    # Output <varname>.clear();
    for vardict in all_vardicts:
        varname = vardict['fullname']
        typename = vardict['branchtype']
        # Clear vectors
        if 'vector' in typename:
            outputline("%s.clear();" % varname)
        else:
            outputline("%s= %s;" % (varname, non_vector_default))

def gen_Branch():
    """Generate Branch() for OutputTree.h"""
    # outputline("#define BRANCH(tree, branch) (tree)->Branch(#branch, &branch);")
    # Output BRANCH(tree, <varname>);
    for vardict in all_vardicts:
        varname = vardict['fullname']
        typename = vardict['branchtype']
        outputline("BRANCH(tree, %s);" % varname)
