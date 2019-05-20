
import MAUS

# Generic Python imports
import sys
import os
import argparse
import math
from math import sqrt
import array

# Third Party library import statements
import json
import event_loader
import analysis
from analysis import tools
from analysis import covariances
from analysis import hit_types
import ROOT

# data cards
#################################
RECON_STATION = 1
RECON_PLANE = 0
SEED_STATION = 1
SEED_PLANE = 0
EXPECTED_STRAIGHT_TRACKPOINTS = 9
EXPECTED_HELIX_TRACKPOINTS = 12
REQUIRE_DATA = True
P_VALUE_CUT = 0.0
MUON_PID = [13, -13]

RECON_TRACKERS = [0, 1]

REQUIRE_ALL_PLANES = True

ANALYSE_REFITS = False

P_MIN = 0.0
P_MAX = 1000.0

PT_MIN = 0.0
PT_MAX = 100.0
PT_BIN = 10
PT_BIN_WIDTH = 10.0

PZ_MIN = 120.0
PZ_MAX = 260.0
PZ_BIN = 14
PZ_BIN_WIDTH = 10.0

X_MIN = -150.0
Y_MIN = -150.0
X_MAX = 150.0
Y_MAX = 150.0

XY_BIN = 300

ALIGNMENT_TOLERANCE = 0.001
RESOLUTION_BINS = 10
EFFICIENCY_BINS = 10

TRACK_ALGORITHM = 1

VIRTUAL_PLANE_DICT = None
INVERSE_PLANE_DICT = {}
TRACKER_PLANE_RADIUS = 100.0

SELECT_EVENTS = False
GOOD_EVENTS = None

SPNPE_CUT = 10.0

######################################################################
def init_plots_data() :
    """
    Setup histogram references & binning
    """

  PZ_BIN = int(((PZ_MAX-PZ_MIN) / PZ_BIN_WIDTH) + 0.5)
  PT_BIN = int(((PT_MAX-PT_MIN) / PT_BIN_WIDTH) + 0.5)
  plot_dict = {'upstream' : {}, 'downstream' : {}}
   
  ####### Upstream & downstream efficiency
  for tracker in [ 'upstream', 'downstream' ] :
    tracker_dict = {}
    tracker_dict['trackpoint_efficiency'] = ROOT.TEfficiency( tracker+'_trackpoint_efficiency', "Track Point Efficiency in P_{z} and P_{#perp}", PZ_BIN, PZ_MIN, PZ_MAX, PT_BIN, PT_MIN, PT_MAX )
    tracker_dict['trackpoint_efficiency_pt'] = ROOT.TEfficiency( tracker+'_trackpoint_efficiency_pt', "Track Point Efficiency in P_{#perp}", PT_BIN, PT_MIN, PT_MAX )
    tracker_dict['trackpoint_efficiency_pz'] = ROOT.TEfficiency( tracker+'_trackpoint_efficiency_pz', "Track Point Efficiency in P_z", PZ_BIN, PZ_MIN, PZ_MAX )

    tracker_dict['spacepoint_efficiency'] = ROOT.TEfficiency( tracker+'_spacepoint_efficiency', "Space Point Efficiency in P_{z} and P_{#perp}", PZ_BIN, PZ_MIN, PZ_MAX, PT_BIN, PT_MIN, PT_MAX )
    tracker_dict['spacepoint_efficiency_pt'] = ROOT.TEfficiency( tracker+'_spacepoint_efficiency_pt', "Space Point Efficiency in P_{#perp}", PT_BIN, PT_MIN, PT_MAX )
    tracker_dict['spacepoint_efficiency_pz'] = ROOT.TEfficiency( tracker+'_spacepoint_efficiency_pz', "Space Point Efficiency in P_z", PZ_BIN, PZ_MIN, PZ_MAX )

    tracker_dict['track_efficiency'] = ROOT.TEfficiency( tracker+'_track_efficiency', "Track Efficiency in P_z and P_{#perp}",  PZ_BIN, PZ_MIN, PZ_MAX, PT_BIN, PT_MIN, PT_MAX )
    tracker_dict['track_efficiency_pt'] = ROOT.TEfficiency(  tracker+'_track_efficiency_pt', "Track Efficiency in P_{#perp}", PT_BIN, PT_MIN, PT_MAX )
    tracker_dict['track_efficiency_pz'] = ROOT.TEfficiency( tracker+'_track_efficiency_pz', "Track Efficiency in P_z", PZ_BIN, PZ_MIN, PZ_MAX )

  ####### Upstream & downstream Purity
    tracker_dict['track_purity'] = ROOT.TEfficiency( tracker+'_track_purity', "Track purity in P_z and P_{#perp}",  PZ_BIN, PZ_MIN, PZ_MAX, PT_BIN, PT_MIN, PT_MAX )
    tracker_dict['track_purity_pt'] = ROOT.TEfficiency(  tracker+'_track_purity_pt', "Track purity in P_{#perp}", PT_BIN, PT_MIN, PT_MAX )
    tracker_dict['track_purity_pz'] = ROOT.TEfficiency( tracker+'_track_purity_pz', "Track purity in P_z", PZ_BIN, PZ_MIN, PZ_MAX )

    tracker_dict['trackpoint_purity'] = ROOT.TEfficiency( tracker+'_trackpoint_purity', "Trackpoint Purity in P_z and P_{#perp}",  PZ_BIN, PZ_MIN, PZ_MAX, PT_BIN, PT_MIN, PT_MAX )
    tracker_dict['trackpoint_purity_pt'] = ROOT.TEfficiency(  tracker+'_trackpoint_purity_pt', "Trackpoint Purity in P_{#perp}", PT_BIN, PT_MIN, PT_MAX )
    tracker_dict['trackpoint_purity_pz'] = ROOT.TEfficiency( tracker+'_trackpoint_purity_pz', "Trackpoint Purity in P_z", PZ_BIN, PZ_MIN, PZ_MAX )

    tracker_dict['spacepoint_purity'] = ROOT.TEfficiency( tracker+'_spacepoint_purity', "Spacepoint Purity in P_z and P_{#perp}",  PZ_BIN, PZ_MIN, PZ_MAX, PT_BIN, PT_MIN, PT_MAX )
    tracker_dict['spacepoint_purity_pt'] = ROOT.TEfficiency(  tracker+'_spacepoint_purity_pt', "Spacepoint Purity in P_{#perp}", PT_BIN, PT_MIN, PT_MAX )
    tracker_dict['spacepoint_purity_pz'] = ROOT.TEfficiency( tracker+'_spacepoint_purity_pz', "Spacepoint Purity in P_z", PZ_BIN, PZ_MIN, PZ_MAX )
  #######  Lost & found plots
    tracker_dict['virtualptpz_found'] = ROOT.TH2F( tracker+'_virtualptpz_found', "Virtuals matched P_z and P_{#perp}",  PZ_BIN, PZ_MIN, PZ_MAX, PT_BIN, PT_MIN, PT_MAX )
    tracker_dict['virtualptpz_missed'] = ROOT.TH2F( tracker+'_virtualptpz_missed', "Virtuals missed P_z and P_{#perp}",  PZ_BIN, PZ_MIN, PZ_MAX, PT_BIN, PT_MIN, PT_MAX )
    tracker_dict['sp_ptpz_found'] = ROOT.TH2F( tracker+'_sp_ptpz_found', "Spacepoint matched P_z and P_{#perp}",  PZ_BIN, PZ_MIN, PZ_MAX, PT_BIN, PT_MIN, PT_MAX )
    tracker_dict['sp_ptpz_missed'] = ROOT.TH2F( tracker+'_sp_ptpz_missed', "Spacepoint missed P_z and P_{#perp}",  PZ_BIN, PZ_MIN, PZ_MAX, PT_BIN, PT_MIN, PT_MAX )
    tracker_dict['sp_xy_found'] = ROOT.TH2F( tracker+'_sp_xy_found', "Spacepoint matched XY",  XY_BIN, X_MIN, X_MAX, XY_BIN, Y_MIN, Y_MAX )
    tracker_dict['sp_xy_missed'] = ROOT.TH2F( tracker+'_sp_xy_missed', "Spacepoint missed XY",  XY_BIN, X_MIN, X_MAX, XY_BIN, Y_MIN, Y_MAX )
    tracker_dict['trackpoint_delta'] = ROOT.TH1F("trackpoint_delta", "trackpoint_delta", 300, 0., 300.)
    tracker_dict['trackpoint_deltaVspt'] = ROOT.TH2F("trackpoint_deltaVspt", "trackpoint_deltaVspt", 100, 0., 100., PT_BIN, PT_MIN, PT_MAX)
    tracker_dict['trackpoint_deltax'] = ROOT.TH1F("trackpoint_deltax", "trackpoint_deltax", 300, -150., 150.)
    tracker_dict['trackpoint_deltay'] = ROOT.TH1F("trackpoint_deltay", "trackpoint_deltay", 300, -150., 150.)
    tracker_dict['trackpoint_deltaz'] = ROOT.TH1F("trackpoint_deltaz", "trackpoint_deltaz", 300, -150., 150.)
    tracker_dict['trackpoint_delta_notfound'] = ROOT.TH1F("trackpoint_delta_notfound", "trackpoint_delta_notfound", 300, 0., 300.)

    tracker_dict['num_tracks'] = ROOT.TH1F("num_tracks", "num_tracks", 10, -0.5, 9.5)

    tracker_dict['novirtual_pt'] = ROOT.TH1F("novirtual_pt", "novirtual_pt", PT_BIN, PT_MIN, PT_MAX)
    tracker_dict['novirtual_pz'] = ROOT.TH1F("novirtual_pz", "novirtual_pz", PZ_BIN, PZ_MIN, PZ_MAX)
    tracker_dict['novirtual_ptpz'] = ROOT.TH2F("novirtual_ptpz", "novirtual_ptpz", PT_BIN, PT_MIN, PT_MAX, PZ_BIN, PZ_MIN, PZ_MAX)
    tracker_dict['novirtual_xy'] = ROOT.TH2F("novirtual_xy", "novirtual_xy", XY_BIN, X_MIN, X_MAX, XY_BIN, Y_MIN, Y_MAX)

    ##### Low-level diagnosis
    plot_dignpe = [[]]
    plot_spnpe = []
    plot_numsp = []
    for station in range (5):
        plot_dignpe.append([])
        name = str(tracker)+"_s"+str(station+1)+"_spnpe"
        plot_spnpe.append(ROOT.TH1F(name, name, 100, 0.0, 100.0))
        tracker_dict[name] = plot_spnpe[station]

        name = str(tracker)+"_s"+str(station+1)+"_numsp"
        plot_numsp.append(ROOT.TH1F(name, name, 10, 0.0, 10.0))
        tracker_dict[name] = plot_numsp[station]
        for plane in range(3):
            name = str(tracker)+"_s"+str(station+1)+"_p"+str(plane)+"_digits"
            plot_dignpe[station].append(ROOT.TH1F(name, name, 301, -0.5, 300.5))
            tracker_dict[name] = plot_dignpe[station][plane]

    plot_dict[tracker] = tracker_dict


  ##### Setup Counters
  data_dict = { 'counters' : {'upstream' : {}, 'downstream' : {} }, \
                                                                  'data' : {} }
  data_dict['counters']['number_events'] = 0

  for tracker in ['upstream', 'downstream'] :
    data_dict['counters'][tracker]['number_virtual'] = 0
    data_dict['counters'][tracker]['missing_virtuals'] = 0

    data_dict['counters'][tracker]['hits'] = 0
    data_dict['counters'][tracker]['hits_within_delta'] = 0

    data_dict['counters'][tracker]['number_tracks'] = 0
    data_dict['counters'][tracker]['number_candidates'] = 0
    data_dict['counters'][tracker]['found_tracks'] = 0
    data_dict['counters'][tracker]['wrong_track_type'] = 0
    data_dict['counters'][tracker]['p_value_cut'] = 0
    data_dict['counters'][tracker]['superfluous_track_events'] = 0
    data_dict['counters'][tracker]['missing_tracks'] = 0
    data_dict['counters'][tracker]['missing_reference_hits'] = 0

    data_dict['counters'][tracker]['momentum_cut'] = 0
    data_dict['counters'][tracker]['gradient_cut'] = 0

    data_dict['counters'][tracker]['found_pairs'] = 0

    data_dict['counters'][tracker]['nvirt'] = 0
    data_dict['counters'][tracker]['recnovirt'] = 0

  return plot_dict, data_dict


######################################################################
def create_virtual_plane_dict(file_reader) :
  """
    Matches up scifitrackpoints to virtual planes to make a lookup dictionary 
  """
  virtual_plane_dict = {}
  for num in range( -15, 0, 1 ) :
    virtual_plane_dict[ num ] = ( -1, (ALIGNMENT_TOLERANCE * 100.0) )
  for num in range( 1, 16, 1 ) :
    virtual_plane_dict[ num ] = ( -1, (ALIGNMENT_TOLERANCE * 100.0) )

  while file_reader.next_event() :
    scifi_event = file_reader.get_event( 'scifi' )
    mc_event = file_reader.get_event( 'mc' )

    tracks = scifi_event.scifitracks()
    for track in tracks :
      if track.tracker() not in RECON_TRACKERS :
        continue
      trackpoints = track.scifitrackpoints()
      for trkpt in trackpoints :
        z_pos = trkpt.pos().z()
        plane_id = analysis.tools.calculate_plane_id(\
                               trkpt.tracker(), trkpt.station(), trkpt.plane())
        #print plane_id,trkpt.tracker(), trkpt.station(), trkpt.plane()

        for vhit_num in xrange(mc_event.GetVirtualHitsSize()) :
          vhit = mc_event.GetAVirtualHit(vhit_num)
          diff = math.fabs(vhit.GetPosition().z() - z_pos)

          if diff < virtual_plane_dict[ plane_id ][1] :
            #print 'adding virtual ',plane_id
            virtual_plane_dict[ plane_id ] = ( vhit.GetStationId(), diff )

    done = True
    for tracker in RECON_TRACKERS :
      for station in [1, 2, 3, 4, 5] :
        for plane in [0, 1, 2] :
          plane_id = analysis.tools.calculate_plane_id( \
                                                      tracker, station, plane )
          if virtual_plane_dict[plane_id][1] > ALIGNMENT_TOLERANCE :
            #print 'TOL: ',plane_id, virtual_plane_dict[plane]
            done = False
    if done :
      break
  else :
    if REQUIRE_ALL_PLANES :
      print
      print virtual_plane_dict
      raise ValueError("Could not locate all virtuals planes")

  file_reader.reset()
  return virtual_plane_dict


######################################################################
def inverse_virtual_plane_dict(virtual_plane_dict) :
  """
    Create the inverse lookup.
  """
  inverse_dict = {}
  for num in range( -15, 0, 1 ) :
    inverse_dict[virtual_plane_dict[num][0]] = num
  for num in range( 1, 16, 1 ) :
    inverse_dict[virtual_plane_dict[num][0]] = num

  return inverse_dict


######################################################################
def get_expected_tracks(mc_event, virtual_plane_dict) :
  upstream_planes = [ virtual_plane_dict[i][0] for i in range(-15, 0)]
  downstream_planes = [ virtual_plane_dict[i][0] for i in range(1, 16)]

  upstream_track = None
  downstream_track = None

  upstream_hits = {}
  downstream_hits = {}

  #print '#virtuals = ',mc_event.GetVirtualHitsSize()
  for vhit_num in xrange(mc_event.GetVirtualHitsSize()) :
    vhit = mc_event.GetAVirtualHit(vhit_num)
    if vhit.GetParticleId() not in MUON_PID :
      continue

    station_id = vhit.GetStationId()
    radius = math.sqrt( vhit.GetPosition().x()**2 + vhit.GetPosition().y()**2 )
    if radius > TRACKER_PLANE_RADIUS : 
      continue

    if station_id in upstream_planes :
      plane_id = INVERSE_PLANE_DICT[station_id]
      upstream_hits[plane_id] = vhit
    if station_id in downstream_planes :
      plane_id = INVERSE_PLANE_DICT[station_id]
      downstream_hits[plane_id] = vhit
  if TRACK_ALGORITHM == 1 :
    if len(upstream_hits) > EXPECTED_HELIX_TRACKPOINTS :
      upstream_track = upstream_hits
    if len(downstream_hits) > EXPECTED_HELIX_TRACKPOINTS :
      downstream_track = downstream_hits
  elif TRACK_ALGORITHM == 0 :
    if len(upstream_hits) > EXPECTED_STRAIGHT_TRACKPOINTS :
      upstream_track = upstream_hits
    if len(downstream_hits) > EXPECTED_STRAIGHT_TRACKPOINTS :
      downstream_track = downstream_hits
  else:
    raise ValueError("Unknown track algorithm found!")

  return upstream_track, downstream_track


######################################################################
def get_found_tracks(scifi_event, plot_dict, data_dict) :
  """
    Find all the single tracks that pass the cuts.
  """
  upstream_tracks = []
  downstream_tracks = []

  tracks = scifi_event.scifitracks()
  nreco_tracks = len(tracks)
  for track in tracks :
    if track.tracker() == 0 :
      tracker = "upstream"
    else :
      tracker = "downstream"

    data_dict['counters'][tracker]['number_tracks'] += 1

    if track.GetAlgorithmUsed() != TRACK_ALGORITHM :
      data_dict['counters'][tracker]['wrong_track_type'] += 1
      continue

    if track.P_value() < P_VALUE_CUT :
      data_dict['counters'][tracker]['p_value_cut'] += 1
      continue

    data_dict['counters'][tracker]['number_candidates'] += 1


    if track.tracker() == 0 :
      upstream_tracks.append(track)
    if track.tracker() == 1 :
      downstream_tracks.append(track)

  if len(upstream_tracks) > 1 :
    data_dict['counters']['upstream']['superfluous_track_events'] += 1
  if len(downstream_tracks) > 1 :
    data_dict['counters']['downstream']['superfluous_track_events'] += 1

  if len(upstream_tracks) == 1 :
    upstream_track = upstream_tracks[0]
    data_dict['counters']['upstream']['found_tracks'] += 1
  else :
    upstream_track = None

  if len(downstream_tracks) == 1 :
    downstream_track = downstream_tracks[0]
    data_dict['counters']['downstream']['found_tracks'] += 1
  else :
    downstream_track = None

  total_tracks_up = len(upstream_tracks)
  total_tracks_down = len(downstream_tracks)
  if total_tracks_up > 0:
    plot_dict['upstream']['num_tracks'].Fill(total_tracks_up)
  if total_tracks_down > 0:
    plot_dict['downstream']['num_tracks'].Fill(total_tracks_down)

  return upstream_track, downstream_track


######################################################################
def make_scifi_mc_pairs(plot_dict, data_dict, virtual_plane_dict, scifi_event, mc_event) :
  """
    Make pairs of SciFiTrackpoints and MC VirtualHits
  """
  paired_hits = []
  paired_seeds = []

  expected_up, expected_down = get_expected_tracks(mc_event, virtual_plane_dict)
  found_up, found_down = get_found_tracks(scifi_event, plot_dict, data_dict)

  downstream_pt = 0.0
  downstream_pz = 0.0

  data_dict['counters']['number_events'] += 1

  for tracker_num, tracker, scifi_track, virtual_track in \
                            [ (0, "upstream", found_up, expected_up), \
                              (1, "downstream", found_down, expected_down) ] :

    if virtual_track is None :
      if scifi_track:
          data_dict['counters'][tracker]['recnovirt'] += 1
          for tp in scifi_track.scifitrackpoints():
              if not tp.has_data():
                  continue
              tt = tp.tracker()
              if tt != tracker_num:
                  continue
              tts = tp.station()
              ttp = tp.plane()
              if tts != 1 or ttp != 0:
                  continue
              rec_px = tp.mom().X()
              rec_py = tp.mom().Y()
              rec_pz = tp.mom().Z()
              rec_pt = math.sqrt(rec_px**2 + rec_py**2)

              rec_x = tp.pos().X()
              rec_y = tp.pos().Y()
              plot_dict[tracker]['track_purity'].Fill(False, rec_pz, rec_pt)
              plot_dict[tracker]['track_purity_pt'].Fill(False, rec_pt)
              plot_dict[tracker]['track_purity_pz'].Fill(False, rec_pt)

              plot_dict[tracker]['novirtual_pt'].Fill(rec_pt)
              plot_dict[tracker]['novirtual_pz'].Fill(rec_pz)
              plot_dict[tracker]['novirtual_ptpz'].Fill(rec_pt,rec_pz)
              plot_dict[tracker]['novirtual_xy'].Fill(rec_x, rec_y)
      continue
    if ANALYSE_REFITS and (scifi_track is None or scifi_track.GetWasRefit() == 0) :
      continue
    data_dict['counters'][tracker]['nvirt'] += 1

    for digit in scifi_event.digits():
        if digit.get_tracker() != tracker_num:
            continue
        s = digit.get_station()-1
        p = digit.get_plane()
        c = digit.get_channel()
        name = str(tracker)+"_s"+str(s+1)+"_p"+str(p)+"_digits"
        plot_dict[tracker][name].Fill(c)

    ref_plane = tools.calculate_plane_id(tracker_num,
                                                    RECON_STATION, RECON_PLANE)
    seed_plane = tools.calculate_plane_id(tracker_num,
                                                    SEED_STATION, SEED_PLANE)
    virtual_pt = 0.0
    virtual_pz = 0.0
    virtual_radius = 0.0
    virtual_hits = 0
    scifi_hits = 0

    seed_virt = None
    reference_virt = None
    reference_scifi = None
    virtual_spacepoints = set()

    for plane in virtual_track :
      if virtual_track[plane] is not None :
        hit = virtual_track[plane]
        px = hit.GetMomentum().x()
        py = hit.GetMomentum().y()
        pt = math.sqrt( px**2 + py**2)
        virtual_pt += pt
        virtual_pz += hit.GetMomentum().z()

        virtual_hits += 1
        if plane == ref_plane :
          reference_virt = virtual_track[plane]
        if plane == seed_plane :
          seed_virt = virtual_track[plane]

        virtual_spacepoints.add(int((abs(plane)-1) / 3))


    #print 'vpt vh: ', virtual_pt, virtual_radius, virtual_hits
    virtual_pt /= virtual_hits
    virtual_pz /= virtual_hits
    virtual_p = math.sqrt( virtual_pt**2 + virtual_pz**2 )
    virtual_radius /= virtual_hits

    virtual_spacepoints = len(virtual_spacepoints)

    if virtual_p > P_MAX or virtual_p < P_MIN :
      data_dict['counters'][tracker]['momentum_cut'] += 1
      continue
    else :
      data_dict['counters'][tracker]['number_virtual'] += 1

    if scifi_track is None :
      plot_dict[tracker]['track_efficiency'].Fill(False, virtual_pz, virtual_pt)
      plot_dict[tracker]['track_efficiency_pt'].Fill(False, virtual_pt)
      plot_dict[tracker]['track_efficiency_pz'].Fill(False, virtual_pz)
      data_dict['counters'][tracker]['missing_tracks'] += 1

      continue # Can't do anything else without a scifi track


    scifi_spacepoints = set()
    for scifi_hit in scifi_track.scifitrackpoints() :
      if scifi_hit.has_data() :
        scifi_hits += 1
        scifi_spacepoints.add( scifi_hit.station()-1 )

        pl_id = analysis.tools.calculate_plane_id(scifi_hit.tracker(), scifi_hit.station(), scifi_hit.plane())

        if scifi_hit.station() == RECON_STATION and scifi_hit.plane() == RECON_PLANE :
          reference_scifi = scifi_hit

    scifi_spacepoints = len(scifi_spacepoints)

    SP_PASSED = True
    nsp = [0]*5
    for sp in scifi_event.spacepoints():
        if sp.get_tracker() != tracker_num:
            continue
        s = sp.get_station()-1
        nsp[s] += 1
        name = str(tracker)+"_s"+str(s+1)+"_spnpe"
        plot_dict[tracker][name].Fill(sp.get_npe())
        if SP_PASSED and sp.get_npe() < SPNPE_CUT:
            SP_PASSED = False
        #### DEBUG
        if sp.get_npe() < 10: 
            print '...... fail ',sp.get_npe()
        print 'pf: ', sp.get_npe(), SP_PASSED
    for s in range(5):
        name = str(tracker)+"_s"+str(s+1)+"_numsp"
        plot_dict[tracker][name].Fill(nsp[s])

    for scifi_hit in scifi_track.scifitrackpoints() :
      if scifi_hit.has_data() :
        scifi_hits += 1

        pl_id = analysis.tools.calculate_plane_id(scifi_hit.tracker(), scifi_hit.station(), scifi_hit.plane())

        if scifi_hit.station() == RECON_STATION and scifi_hit.plane() == RECON_PLANE :
          reference_scifi = scifi_hit

        if scifi_hit.tracker() != tracker_num:
            continue
        if scifi_hit.station() != RECON_STATION or scifi_hit.plane() != RECON_PLANE:
            continue
        rec_px = scifi_hit.mom().X()
        rec_py = scifi_hit.mom().Y()
        rec_pz = scifi_hit.mom().Z()
        rec_pt = math.sqrt(rec_px**2 + rec_py**2)
        plot_dict[tracker]['track_purity'].Fill(True, rec_pz, rec_pt)
        plot_dict[tracker]['track_purity_pt'].Fill(True, rec_pt)
        plot_dict[tracker]['track_purity_pz'].Fill(True, rec_pt)

    if not SP_PASSED:
        continue

    plot_dict[tracker]['track_efficiency'].Fill(True, virtual_pz, virtual_pt)
    plot_dict[tracker]['track_efficiency_pt'].Fill(True, virtual_pt)
    plot_dict[tracker]['track_efficiency_pz'].Fill(True, virtual_pz)

    if scifi_hits >= virtual_hits :
      for i in range(virtual_hits) :
        plot_dict[tracker]['trackpoint_efficiency'].Fill(True, virtual_pz, virtual_pt)
        plot_dict[tracker]['trackpoint_efficiency_pt'].Fill(True, virtual_pt)
        plot_dict[tracker]['trackpoint_efficiency_pz'].Fill(True, virtual_pz)
    else :
      for i in range( virtual_hits - scifi_hits ) :
        plot_dict[tracker]['trackpoint_efficiency'].Fill(False, virtual_pz, virtual_pt)
        plot_dict[tracker]['trackpoint_efficiency_pt'].Fill(False, virtual_pt)
        plot_dict[tracker]['trackpoint_efficiency_pz'].Fill(False, virtual_pz)
      for i in range( scifi_hits ) :
        plot_dict[tracker]['trackpoint_efficiency'].Fill(True, virtual_pz, virtual_pt)
        plot_dict[tracker]['trackpoint_efficiency_pt'].Fill(True, virtual_pt)
        plot_dict[tracker]['trackpoint_efficiency_pz'].Fill(True, virtual_pz)

    if scifi_spacepoints >= virtual_spacepoints :
      for i in range(virtual_hits) :
        plot_dict[tracker]['spacepoint_efficiency'].Fill(True, virtual_pz, virtual_pt)
        plot_dict[tracker]['spacepoint_efficiency_pt'].Fill(True, virtual_pt)
        plot_dict[tracker]['spacepoint_efficiency_pz'].Fill(True, virtual_pz)
    else :
      for i in range( virtual_spacepoints - scifi_spacepoints ) :
        plot_dict[tracker]['spacepoint_efficiency'].Fill(False, virtual_pz, virtual_pt)
        plot_dict[tracker]['spacepoint_efficiency_pt'].Fill(False, virtual_pt)
        plot_dict[tracker]['spacepoint_efficiency_pz'].Fill(False, virtual_pz)
      for i in range( scifi_spacepoints ) :
        plot_dict[tracker]['spacepoint_efficiency'].Fill(True, virtual_pz, virtual_pt)
        plot_dict[tracker]['spacepoint_efficiency_pt'].Fill(True, virtual_pt)
        plot_dict[tracker]['spacepoint_efficiency_pz'].Fill(True, virtual_pz)

    spacepoint_set = set()
    for scifi_hit in scifi_track.scifitrackpoints() :
      if scifi_hit.has_data() :
        scifi_plane = analysis.tools.calculate_plane_id(scifi_hit.tracker(), scifi_hit.station(), scifi_hit.plane())

        if scifi_plane in virtual_track :
          virtual_hit = virtual_track[scifi_plane]
          vx = virtual_hit.GetPosition().x()
          vy = virtual_hit.GetPosition().y()
          vz = virtual_hit.GetPosition().z()
          sx = scifi_hit.pos().x()
          sy = scifi_hit.pos().y()
          sz = scifi_hit.pos().z()
          delta = math.sqrt((sx-vx)**2 + (sy-vy)**2)
          plot_dict[tracker]['trackpoint_deltax'].Fill(sx-vx)
          plot_dict[tracker]['trackpoint_deltay'].Fill(sy-vy)
          plot_dict[tracker]['trackpoint_deltaz'].Fill(sz-vz)
          plot_dict[tracker]['trackpoint_delta'].Fill(delta)
          plot_dict[tracker]['trackpoint_deltaVspt'].Fill(delta, virtual_pt)
         
          if scifi_hit.station() == RECON_STATION and scifi_hit.plane() == RECON_PLANE:
           data_dict['counters'][tracker]['hits'] += 1
           if delta < 5.0 :
            data_dict['counters'][tracker]['hits_within_delta'] += 1
            plot_dict[tracker]["trackpoint_purity"].Fill(True, virtual_pz, virtual_pt)
            plot_dict[tracker]["trackpoint_purity_pt"].Fill(True, virtual_pt)
            plot_dict[tracker]["trackpoint_purity_pz"].Fill(True, virtual_pz)
           else :
            plot_dict[tracker]['trackpoint_delta_notfound'].Fill(delta)
            plot_dict[tracker]["trackpoint_purity"].Fill(False, virtual_pz, virtual_pt)
            plot_dict[tracker]["trackpoint_purity_pt"].Fill(False, virtual_pt)
            plot_dict[tracker]["trackpoint_purity_pz"].Fill(False, virtual_pz)

          if not (scifi_hit.station() in spacepoint_set) :
            spacepoint_set.add(scifi_hit.station())
            if delta < 5.0 :
              plot_dict[tracker]["spacepoint_purity"].Fill(True, virtual_pz, virtual_pt)
              plot_dict[tracker]["spacepoint_purity_pt"].Fill(True, virtual_pt)
              plot_dict[tracker]["spacepoint_purity_pz"].Fill(True, virtual_pz)
            else :
              plot_dict[tracker]["spacepoint_purity"].Fill(False, virtual_pz, virtual_pt)
              plot_dict[tracker]["spacepoint_purity_pt"].Fill(False, virtual_pt)
              plot_dict[tracker]["spacepoint_purity_pz"].Fill(False, virtual_pz)


    if reference_virt is None :
      data_dict['counters'][tracker]['missing_virtuals'] += 1

    if reference_scifi is None :
      data_dict['counters'][tracker]['missing_reference_hits'] += 1

    if reference_virt is not None and reference_scifi is not None :
      paired_hits.append( (reference_scifi, reference_virt) )
      data_dict['counters'][tracker]['found_pairs'] += 1

    if seed_virt is not None and scifi_track is not None :
      paired_seeds.append( (scifi_track, seed_virt))

  return paired_hits, paired_seeds


######################################################################
def analyse_plots(plot_dict, data_dict) :

  for tracker in [ 'upstream', 'downstream' ] :
    plot_dict[tracker]['trackpoint_efficiency'] = plot_dict[tracker]['trackpoint_efficiency'].CreateHistogram()
    plot_dict[tracker]['trackpoint_efficiency_pt'] = plot_dict[tracker]['trackpoint_efficiency_pt'].CreateGraph()
    plot_dict[tracker]['trackpoint_efficiency_pz'] = plot_dict[tracker]['trackpoint_efficiency_pz'].CreateGraph()

    plot_dict[tracker]['spacepoint_efficiency'] = plot_dict[tracker]['spacepoint_efficiency'].CreateHistogram()
    plot_dict[tracker]['spacepoint_efficiency_pt'] = plot_dict[tracker]['spacepoint_efficiency_pt'].CreateGraph()
    plot_dict[tracker]['spacepoint_efficiency_pz'] = plot_dict[tracker]['spacepoint_efficiency_pz'].CreateGraph()

    plot_dict[tracker]['track_efficiency'] = plot_dict[tracker]['track_efficiency'].CreateHistogram()
    plot_dict[tracker]['track_efficiency_pt'] = plot_dict[tracker]['track_efficiency_pt'].CreateGraph()
    plot_dict[tracker]['track_efficiency_pz'] = plot_dict[tracker]['track_efficiency_pz'].CreateGraph()
                                                                                                   
    plot_dict[tracker]['track_purity'] = plot_dict[tracker]['track_purity'].CreateHistogram()
    plot_dict[tracker]['track_purity_pt'] = plot_dict[tracker]['track_purity_pt'].CreateGraph()
    plot_dict[tracker]['track_purity_pz'] = plot_dict[tracker]['track_purity_pz'].CreateGraph()
                                                                                                   
    plot_dict[tracker]['trackpoint_purity'] = plot_dict[tracker]['trackpoint_purity'] .CreateHistogram()
    plot_dict[tracker]['trackpoint_purity_pt'] = plot_dict[tracker]['trackpoint_purity_pt'].CreateGraph()
    plot_dict[tracker]['trackpoint_purity_pz'] = plot_dict[tracker]['trackpoint_purity_pz'].CreateGraph()

    plot_dict[tracker]['spacepoint_purity'] = plot_dict[tracker]['spacepoint_purity'] .CreateHistogram()
    plot_dict[tracker]['spacepoint_purity_pt'] = plot_dict[tracker]['spacepoint_purity_pt'].CreateGraph()
    plot_dict[tracker]['spacepoint_purity_pz'] = plot_dict[tracker]['spacepoint_purity_pz'].CreateGraph()

    print data_dict['counters'][tracker]['hits'], data_dict['counters'][tracker]['hits_within_delta']


######################################################################
if __name__ == "__main__" : 
  ROOT.gROOT.SetBatch( True )
  ROOT.gErrorIgnoreLevel = ROOT.kError

  parser = argparse.ArgumentParser( description='An example script showing '+\
      'some basic data extraction and analysis routines' )

  parser.add_argument( 'maus_root_files', nargs='+', help='List of MAUS '+\
                  'output root files containing reconstructed straight tracks')

  parser.add_argument( '-N', '--max_num_events', type=int, \
                                   help='Maximum number of events to analyse.')

  parser.add_argument( '-O', '--output_filename', \
            default='tracker_purity_plots', help='Set the output filename')

  parser.add_argument( '-D', '--output_directory', \
                                 default='./', help='Set the output directory')

  parser.add_argument( '-V', '--virtual_plane_dictionary', default=None, \
                   help='Specify a json file containing a dictionary of the '+\
                                                       'virtual plane lookup' )

  parser.add_argument( '-P', '--print_plots', action='store_true', \
                        help="Flag to save the plots as individual pdf files" )

  parser.add_argument( '--cut_p_value', type=float, default=0.0, \
  help="Specify the P-Value below which tracks are removed from the analysis" )

  parser.add_argument( '--track_algorithm', type=int, default=1, \
                          help="Specify the track reconstruction algorithm. "+\
                             "1 for Helical Tracks and 0 for Straight Tracks" )

  parser.add_argument( '--pz_bin', type=float, default=PZ_BIN_WIDTH, \
             help="Specify the size of the Pz bins which are used to select "+\
                     "particles for the reconstruction of optical functions." )
  parser.add_argument( '--pz_window', type=float, nargs=2, \
        default=[PZ_MIN, PZ_MAX], help="Specify the range of Pz to consider "+\
                               "for the reconstruction of optical functions." )

  parser.add_argument( '--pt_bin', type=float, default=PT_BIN_WIDTH, \
             help="Specify the size of the Pt bins which are used to select "+\
                     "particles for the reconstruction of optical functions." )
  parser.add_argument( '--pt_window', type=float, nargs=2, \
        default=[PT_MIN, PT_MAX], help="Specify the range of Pt to consider "+\
                               "for the reconstruction of optical functions." )

  parser.add_argument( '--trackers', type=int, default=RECON_TRACKERS, \
                          nargs='+', help="Specifies the trackers to analyse" )

  parser.add_argument( '--p_window', type=float, nargs=2, \
             default=[P_MIN, P_MAX], help="Specify the range of the total " + \
                                         "momentum to consider for analysis." )

  parser.add_argument( '--selection_file', default=None, \
                 help='Name of a JSON file containing the events to analyses' )


  parser.add_argument( '--analyse_refits', action='store_true', \
                       help='Only reconstruct the tracks flagged as refitted' )


  try :
    namespace = parser.parse_args()

    P_VALUE_CUT = namespace.cut_p_value
    TRACK_ALGORITHM = namespace.track_algorithm

    RECON_TRACKERS = namespace.trackers

    ANALYSE_REFITS = namespace.analyse_refits

    P_MIN = namespace.p_window[0]
    P_MAX = namespace.p_window[1]

    PZ_MIN = namespace.pz_window[0]
    PZ_MAX = namespace.pz_window[1]
    PZ_BIN_WIDTH = namespace.pz_bin

    PT_MIN = namespace.pt_window[0]
    PT_MAX = namespace.pt_window[1]
    PT_BIN_WIDTH = namespace.pt_bin

    if namespace.selection_file is not None :
      SELECT_EVENTS = True
      with open(namespace.selection_file, 'r') as infile :
        GOOD_EVENTS = json.load(infile)
    else :
      SELECT_EVENTS = False

    if namespace.virtual_plane_dictionary is not None :
      VIRTUAL_PLANE_DICT = analysis.tools.load_virtual_plane_dict( \
                                           namespace.virtual_plane_dictionary )
  except BaseException as ex:
    raise
  else :
##### 1. Load MAUS globals and geometry. - NOT NECESSARY AT PRESENT
    # geom = load_tracker_geometry(namespace.configuration_file)

##### 2. Intialise plots ######################################################
    print
    sys.stdout.write( "\n- Initialising Plots : Running\r" )
    sys.stdout.flush()
    plot_dict, data_dict = init_plots_data()
    sys.stdout.write(   "- Initialising Plots : Done   \n" )

    file_reader = event_loader.maus_reader(namespace.maus_root_files)

##### 3. Initialise Plane Dictionary ##########################################
    if VIRTUAL_PLANE_DICT is None :
      sys.stdout.write( "\n- Finding Virtual Planes : Running\r" )
      sys.stdout.flush()
      virtual_plane_dictionary = create_virtual_plane_dict(file_reader)
      VIRTUAL_PLANE_DICT = virtual_plane_dictionary
      sys.stdout.write(   "- Finding Virtual Planes : Done   \n" )

    INVERSE_PLANE_DICT = inverse_virtual_plane_dict(VIRTUAL_PLANE_DICT)
    file_reader.select_events(GOOD_EVENTS)
    file_reader.set_max_num_events(namespace.max_num_events)
    file_reader.set_print_progress('spill')

##### 4. Load Events ##########################################################
    print "\n- Loading Spills...\n"
    try :
      while file_reader.next_selected_event() :

        try :
          scifi_event = file_reader.get_event( 'scifi' )
          mc_event = file_reader.get_event( 'mc' )
          tof_event = file_reader.get_event( 'tof' )
          tofsps = tof_event.GetTOFEventSpacePoint()
          tof1_sps = tofsps.GetTOF1SpacePointArray()
          if len(tof1_sps) != 1:
              continue


##### 5. Extract tracks and Fill Plots ########################################

          paired_hits, seed_pairs = make_scifi_mc_pairs(plot_dict, data_dict, \
                                     VIRTUAL_PLANE_DICT, scifi_event, mc_event)
        except ValueError as ex :
          print "An Error Occured: " + str(ex)
          print "Skipping Event: " +\
                str(file_reader.get_current_event_number()) + " In Spill: " + \
                str(file_reader.get_current_spill_number()) + " In File: " + \
                str(file_reader.get_current_filenumber()) + "\n"
          continue

    except KeyboardInterrupt :
      print
      print " ###  Keyboard Interrupt  ###"
      print
    print "- {0:0.0f} Spills Loaded                                 ".format( \
                                            file_reader.get_total_num_spills())
##### 6. Analysing Plots ######################################################
    print"\n- Analysing Data...\n"

    analyse_plots(plot_dict, data_dict)

##### 7. Saving Plots and Data ################################################

    sys.stdout.write( "\n- Saving Plots and Data : Running\r" )
    sys.stdout.flush()
    save_pretty(plot_dict, namespace.output_directory )

    filename = os.path.join(namespace.output_directory, \
                                                     namespace.output_filename)
    analysis.tools.save_plots(plot_dict, filename+'.root')
    sys.stdout.write(   "- Saving Plots and Data : Done   \n" )

  print 
  print "Complete."
  print

