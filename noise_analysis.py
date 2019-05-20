import MAUS

import sys
import os
import argparse

import math
import json
import event_loader
import analysis
from analysis import tools
from analysis import covariances
from analysis import hit_types
import ROOT

###############################################
# Cards, counters, flags
#########################
NPE_CUT = 2.0

NEIGHBOUR_PT_CUT = 2000.0

TOF01_CUT_LOW = 28.0
TOF01_CUT_HIGH = 31.0
TOF12_CUT_LOW = 32.0
TOF12_CUT_HIGH = 40.0

TOF01_PLOT = ROOT.TH1F('tof01_plot', '', 200, 20.0, 40.0 )
TOF12_PLOT = ROOT.TH1F('tof12_plot', '', 400, 20.0, 60.0 )

TOF0_PASSED = False
TOF1_PASSED = False
TOF2_PASSED = False
TOF12_PASSED = False
TOF01_PASSED = False

#################################
def passes_cuts( tof_event ):
  """
  Check if a TOF event passes the selection requirements
  Checks # hits, # reconstructed space-points, reconstructed flight times
  """

  global TOF0_PASSED
  global TOF1_PASSED
  global TOF2_PASSED
  global TOF12_PASSED
  global TOF01_PASSED

  TOF0_PASSED = False
  TOF1_PASSED = False
  TOF2_PASSED = False
  TOF12_PASSED = False
  TOF01_PASSED = False

  event_spacepoints = tof_event.GetTOFEventSpacePoint()

  tof0_sp_size = event_spacepoints.GetTOF0SpacePointArraySize()
  tof1_sp_size = event_spacepoints.GetTOF1SpacePointArraySize()
  tof2_sp_size = event_spacepoints.GetTOF2SpacePointArraySize()

  if tof0_sp_size == 1 :
    TOF0_PASSED = True

  if tof1_sp_size == 1:
    TOF1_PASSED = True
  
  if tof2_sp_size == 1 :
    TOF2_PASSED = True

  tof0_time = -999
  tof1_time = -999
  tof2_time = -999

  if TOF0_PASSED:
    tof0_sp = event_spacepoints.GetTOF0SpacePointArrayElement(0)
    tof0_time = tof0_sp.GetTime()
  if TOF1_PASSED:
    tof1_sp = event_spacepoints.GetTOF1SpacePointArrayElement(0)
    tof1_time = tof1_sp.GetTime()
  if TOF2_PASSED:
    tof2_sp = event_spacepoints.GetTOF2SpacePointArrayElement(0)
    tof2_time = tof2_sp.GetTime()


  diff_0_1 = 999
  diff_1_2 = 999

  if tof0_time > -999 and tof1_time > -999:
    diff_0_1 = tof1_sp.GetTime() - tof0_sp.GetTime()
  if tof1_time > -999 and tof2_time > -999:
    diff_1_2 = tof2_sp.GetTime() - tof1_sp.GetTime()

  TOF01_PLOT.Fill(diff_0_1)
  TOF12_PLOT.Fill(diff_1_2)

  if diff_0_1 >= TOF01_CUT_LOW and diff_0_1 <= TOF01_CUT_HIGH :
    TOF01_PASSED = True
  if diff_1_2 >= TOF12_CUT_LOW and diff_1_2 <= TOF12_CUT_HIGH :
    TOF12_PASSED = True

  return True


###############################################
class Tracker_Analyser(object) :
    """
    Primary class
    Handles initialization, lookup, analysis, plotting
    """

  ###############################################
  def __init__(self, tracker_no) :
    self.__tracker = tracker_no

    self.__total_event_counter = 0
    self.__event_counter = 0
    self.__good_event_counter = 0
    self.__bad_event_counter = 0
    self.__digit_counter = 0
    self.__cluster_counter = 0
    self.__spacepoint_counter = 0
    self.__spacepoint_CUT_counter = 0
    self.__helical_counter = 0
    self.__straight_counter = 0
    self.__helical_events = 0
    self.__straight_events = 0

    self.__expected_track_events = 0
    self.__no_expected_track_events = 0
    self.__expected_missing_track_events = 0

    self.__new_select_missed = 0
    self.__old_select_missed = 0

    self.__new_missed = 0
    self.__old_missed = 0
    self.__new_found = 0
    self.__old_found = 0

    self.__neighbour_found = 0
    self.__neighbour_missed = 0
    self.__neighbour_expected_missed = 0

    self.__num_expected_other_tracker = 0
    self.__no_expected_other = 0
    self.__missing_here = 0

    self.__plot_clnpe = [[]]
    self.__plot_clnpe_expected = [[]]
    self.__plot_clnpe_no_expected = [[]]
    self.__plot_dignpe = [[]]
    self.__plot_dignpe_expected = [[]]
    self.__plot_dignpeVchan = [[]]
    self.__plot_digadc = [[]]
    self.__plot_digits = [[]]
    self.__plot_spnpe_used = []
    self.__plot_spnpe_notused = []
    self.__plot_spxy_used = []
    self.__plot_spxy_notused = []
    self.__plot_spxy_cutevt = []
    self.__plot_spnpe_cutevt = []
    for station in range (5):
        name = str(self.__tracker)+"_s"+str(station+1)+"_spnpe_used"
        self.__plot_spnpe_used.append(ROOT.TH1F(name, name, 100, 0., 100.))

        name = str(self.__tracker)+"_s"+str(station+1)+"_spnpe_notused"
        self.__plot_spnpe_notused.append(ROOT.TH1F(name, name, 100, 0., 100.))

        name = str(self.__tracker)+"_s"+str(station+1)+"_spxy_used"
        self.__plot_spxy_used.append(ROOT.TH2F(name, name, 100, -150., 150., 100, -150., 150.))

        name = str(self.__tracker)+"_s"+str(station+1)+"_spxy_notused"
        self.__plot_spxy_notused.append(ROOT.TH2F(name, name, 100, -150., 150., 100, -150., 150.))

        name = str(self.__tracker)+"_s"+str(station+1)+"_spxy_cutevt"
        self.__plot_spxy_cutevt.append(ROOT.TH2F(name, name, 100, -150., 150., 100, -150., 150.))

        name = str(self.__tracker)+"_s"+str(station+1)+"_spnpe_cutevt"
        self.__plot_spnpe_cutevt.append(ROOT.TH1F(name, name, 100, 0., 100.))

        self.__plot_clnpe.append([])
        self.__plot_clnpe_expected.append([])
        self.__plot_clnpe_no_expected.append([])
        self.__plot_dignpe.append([])
        self.__plot_dignpe_expected.append([])
        self.__plot_dignpeVchan.append([])
        self.__plot_digadc.append([])
        self.__plot_digits.append([])
        for plane in range(3):
            name = str(self.__tracker)+"_s"+str(station+1)+"_p"+str(plane)+"_cluster_npe"
            self.__plot_clnpe[station].append(ROOT.TH1F(name, name, 100, 0., 100.))
            name = str(self.__tracker)+"_s"+str(station+1)+"_p"+str(plane)+"_cluster_npe_expected"
            self.__plot_clnpe_expected[station].append(ROOT.TH1F(name, name, 100, 0., 100.))
            name = str(self.__tracker)+"_s"+str(station+1)+"_p"+str(plane)+"_cluster_npe_not_expected"
            self.__plot_clnpe_no_expected[station].append(ROOT.TH1F(name, name, 100, 0., 100.))

            name = str(self.__tracker)+"_s"+str(station+1)+"_p"+str(plane)+"_dig_npe"
            self.__plot_dignpe[station].append(ROOT.TH1F(name, name, 100, 0., 100.))
            name = str(self.__tracker)+"_s"+str(station+1)+"_p"+str(plane)+"_dig_npe_expected"
            self.__plot_dignpe_expected[station].append(ROOT.TH1F(name, name, 100, 0., 100.))

            name = str(self.__tracker)+"_s"+str(station+1)+"_p"+str(plane)+"_dig_npeVchan"
            self.__plot_dignpeVchan[station].append(ROOT.TH2F(name, name, 301, -0.5, 300.5, 100, 0., 100.))

            name = str(self.__tracker)+"_s"+str(station+1)+"_p"+str(plane)+"_dig_adc"
            self.__plot_digadc[station].append(ROOT.TH1F(name, name, 100, 0., 4100.))

            name = str(self.__tracker)+"_s"+str(station+1)+"_p"+str(plane)+"_digits"
            self.__plot_digits[station].append(ROOT.TH1F(name, name, 301, -0.5, 300.5))


    self.__plot_cluster_npe = ROOT.TH1F(str(self.__tracker)+"_cluster_npe", "", 100, 0.0, 100.0)
    self.__plot_cluster_npe_expected = ROOT.TH1F(str(self.__tracker)+"_cluster_npe_expected", "", 100, 0.0, 100.0)
    self.__plot_cluster_npe_notexpected = ROOT.TH1F(str(self.__tracker)+"_cluster_npe_notexpected", "", 100, 0.0, 100.0)
    self.__plot_cluster_npe_neighbours = ROOT.TH1F(str(self.__tracker)+"_cluster_npe_neighbours", "", 100, 0.0, 100.0)
    self.__plot_cluster_npe_noneighbours = ROOT.TH1F(str(self.__tracker)+"_cluster_npe_noneighbours", "", 100, 0.0, 100.0)

    self.__plot_spacepoint_npe = ROOT.TH1F(str(self.__tracker)+"_spacepoint_npe", "", 100, 0.0, 100.0)
    self.__plot_cutevt_spacepoint_npe = ROOT.TH1F(str(self.__tracker)+"_cutevt_spacepoint_npe", "", 100, 0.0, 100.0)
    self.__plot_spacepoints = ROOT.TH1F(str(self.__tracker)+"_spacepoint", "", 101, -0.5, 100.5)
    self.__plot_spacepoints_CUT = ROOT.TH1F(str(self.__tracker)+"_spacepoint_CUT", "", 101, -0.5, 100.5)

    self.__plot_found_spacepoint_npe = ROOT.TH1F(str(self.__tracker)+"_found_spacepoint_npe", "", 100, 0.0, 100.0)
    self.__plot_unused_spacepoint_npe = ROOT.TH1F(str(self.__tracker)+"_unused_spacepoint_npe", "", 100, 0.0, 100.0)

    self.__plot_expected_spnpe = ROOT.TH1F(str(self.__tracker)+"_expected_spnpe", "", 100, 0.0, 100.0)
    self.__plot_notexpected_spnpe = ROOT.TH1F(str(self.__tracker)+"_notexpected_spnpe", "", 100, 0.0, 100.0)
    self.__plot_neighbours_spnpe = ROOT.TH1F(str(self.__tracker)+"_neighbours_spnpe", "", 100, 0.0, 100.0)
    self.__plot_noneighbours_spnpe = ROOT.TH1F(str(self.__tracker)+"_noneighbours_spnpe", "", 100, 0.0, 100.0)

    self.__plot_cutevt_expected_spnpe = ROOT.TH1F(str(self.__tracker)+"_cutevt_expected_spnpe", "", 100, 0.0, 100.0)
    self.__plot_cutevt_notexpected_spnpe = ROOT.TH1F(str(self.__tracker)+"_cutevt_notexpected_spnpe", "", 100, 0.0, 100.0)
    self.__plot_cutevt_neighbours_spnpe = ROOT.TH1F(str(self.__tracker)+"_cutevt_neighbours_spnpe", "", 100, 0.0, 100.0)
    self.__plot_cutevt_noneighbours_spnpe = ROOT.TH1F(str(self.__tracker)+"_cutevt_noneighbours_spnpe", "", 100, 0.0, 100.0)

    self.__plot_expected_spxy = ROOT.TH2F(str(self.__tracker)+"_expected_spxy", "", 100, -150.0, 150.0, 100, -150.0, 150.)
    self.__plot_notexpected_spxy = ROOT.TH2F(str(self.__tracker)+"_notexpected_spxy", "", 100, 150.0, 150.0, 100, -150.0, 150.)
    self.__plot_neighbours_spxy = ROOT.TH2F(str(self.__tracker)+"_neighbours_spxy", "", 100, -150.0, 150.0, 100, -150.0, 150.)
    self.__plot_noneighbours_spxy = ROOT.TH2F(str(self.__tracker)+"_noneighbours_spxy", "", 100, 150.0, 150.0, 100, -150.0, 150.)

    self.__plot_cutevt_expected_spxy = ROOT.TH2F(str(self.__tracker)+"_cutevt_expected_spxy", "", 100, -150.0, 150.0, 100, -150.0, 150.)
    self.__plot_cutevt_notexpected_spxy = ROOT.TH2F(str(self.__tracker)+"_cutevt_notexpected_spxy", "", 100, 150.0, 150.0, 100, -150.0, 150.)
    self.__plot_cutevt_neighbours_spxy = ROOT.TH2F(str(self.__tracker)+"_cutevt_neighbours_spxy", "", 100, -150.0, 150.0, 100, -150.0, 150.)
    self.__plot_cutevt_noneighbours_spxy = ROOT.TH2F(str(self.__tracker)+"_cutevt_noneighbours_spxy", "", 100, 150.0, 150.0, 100, -150.0, 150.)

    self.__plot_spacepoints_CUT_new = ROOT.TH1F(str(self.__tracker)+"_spacepoint_CUT_new", "", 101, -0.5, 100.5)
    self.__plot_spacepoints_CUT_old = ROOT.TH1F(str(self.__tracker)+"_spacepoint_CUT_old", "", 101, -0.5, 100.5)

    self.__plot_spacepoints_station_good = ROOT.TH2F(str(self.__tracker)+"_spacepoint_station_good", "", 5, 0.5, 5.5, 10, -0.5, 9.5)
    self.__plot_spacepoints_station_old = ROOT.TH2F(str(self.__tracker)+"_spacepoint_station_old", "", 5, 0.5, 5.5, 10, -0.5, 9.5)
    self.__plot_spacepoints_station_new = ROOT.TH2F(str(self.__tracker)+"_spacepoint_station_new", "", 5, 0.5, 5.5, 10, -0.5, 9.5)

    self.__plot_chisq_ndf = ROOT.TH1F(str(self.__tracker)+"_chisq_ndf", "", 500, 0.0, 100.0 )

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

    PZ_BIN = int(((PZ_MAX-PZ_MIN) / PZ_BIN_WIDTH) + 0.5)
    PT_BIN = int(((PT_MAX-PT_MIN) / PT_BIN_WIDTH) + 0.5)
    self.__plot_track_efficiency = ROOT.TEfficiency( str(self.__tracker)+'_track_efficiency', "Track Efficiency in P_z and P_{#perp}",  PZ_BIN, PZ_MIN, PZ_MAX, PT_BIN, PT_MIN, PT_MAX )
    self.__plot_track_efficiency_pt = ROOT.TEfficiency(  str(self.__tracker)+'_track_efficiency_pt', "Track Efficiency in P_{#perp}", PT_BIN, PT_MIN, PT_MAX )
    self.__plot_track_efficiency_pz = ROOT.TEfficiency( str(self.__tracker)+'_track_efficiency_pz', "Track Efficiency in P_z", PZ_BIN, PZ_MIN, PZ_MAX )



  ###############################################
  def analyse_event(self, scifi_event, passed) :
    """
    INPUT: reference to a scifi event 
    OUTPUT: pass/fail 
    """
    digits = scifi_event.digits()
    clusters = scifi_event.clusters()
    spacepoints = scifi_event.spacepoints()

    helical_tracks = scifi_event.helicalprtracks()
    straight_tracks = scifi_event.straightprtracks()

    kalman_tracks = scifi_event.scifitracks()

    has_neighbours = False

    tracks_expected_other = False
    if self.__tracker == 0 :
      tracks_expected = scifi_event.get_expected_track_upstream()
      tracks_expected_other = scifi_event.get_expected_track_downstream()
    else :
      tracks_expected = scifi_event.get_expected_track_downstream()
      tracks_expected_other = scifi_event.get_expected_track_upstream()

    num_digits = 0
    num_clusters = 0
    num_spacepoints = 0
    num_spacepoints_CUT = 0
    num_helicals = 0
    num_straights = 0

    num_other_helicals = 0
    num_other_straights = 0
    num_other_spacepoints = 0
    num_other_spacepoints_CUT = 0

    num_neighbours = 0
    spacepoints_station = { 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0 }

    self.__total_event_counter += 1

    ###### helicals before any cuts
    for helical in helical_tracks :
      if helical.get_tracker() == self.__tracker :
        num_helicals += 1
      elif helical.get_reference_momentum().Pt() < NEIGHBOUR_PT_CUT :
        has_neighbours = True
        num_neighbours += 1
      if helical.get_tracker() != self.__tracker:
          num_other_helicals += 1

    ###### straights before any cuts
    for straight in straight_tracks :
      if straight.get_tracker() == self.__tracker :
        num_straights += 1
      else:
        num_other_straights += 1

    # total tracks before any cuts
    total_tracks = num_straights + num_helicals
    other_total_tracks = num_other_helicals + num_other_straights

    # digits before any cuts
    for digit in digits :
      if digit.get_tracker() == self.__tracker :
        num_digits += 1

    # clusters before any cuts
    for cluster in clusters :
      if cluster.get_tracker() == self.__tracker :
        num_clusters += 1

    # spacepoints before any cuts
    for spacepoint in spacepoints :
      if spacepoint.get_tracker() == self.__tracker :
        num_spacepoints += 1
        spacepoints_station[ spacepoint.get_station() ] += 1.0
        if spacepoint.get_npe() > NPE_CUT :
          num_spacepoints_CUT += 1
      else: # other tracker
        num_other_spacepoints += 1
        if spacepoint.get_npe() > NPE_CUT :
          num_other_spacepoints_CUT += 1

    if num_spacepoints_CUT > 3 :
      self.__plot_spacepoints_CUT_new.Fill( num_spacepoints_CUT )
      for key, value in spacepoints_station.iteritems() :
        self.__plot_spacepoints_station_new.Fill( key, value )
    if total_tracks > 0 :
      self.__plot_spacepoints_CUT_old.Fill( num_spacepoints_CUT )
      for key, value in spacepoints_station.iteritems() :
        self.__plot_spacepoints_station_old.Fill( key, value )


    if self.__tracker == 0:
        # upstream, check if TOF1 & downstream tracker
        EVENT_PASSED = TOF0_PASSED and TOF1_PASSED and TOF01_PASSED
    else:
        EVENT_PASSED = TOF0_PASSED and TOF1_PASSED and TOF01_PASSED and TOF2_PASSED and TOF12_PASSED and other_total_tracks==1

    if EVENT_PASSED:
      ###### DEBUG PRINT
      # print 'tracks: other: passed',total_tracks, other_total_tracks, num_neighbours, EVENT_PASSED, num_spacepoints, num_spacepoints_CUT, num_other_spacepoints, num_other_spacepoints_CUT
      self.__event_counter += 1

      rad = 0.0
      RAD_CUT = False
      for track in kalman_tracks :
        tps = track.scifitrackpoints()
        RAD_CUT = False
        for tp in tps:
            if RAD_CUT:
                continue
            x = tp.pos().X()
            y = tp.pos().Y()
            rad = math.sqrt(x**2 + y**2)
            if rad > 100:
                #print '> ',rad
                RAD_CUT = True

      if RAD_CUT:
          return
      self.__helical_counter += num_helicals
      self.__straight_counter += num_straights
      self.__digit_counter += num_digits
      self.__cluster_counter += num_clusters
      self.__spacepoint_counter += num_spacepoints

      if has_neighbours :
        self.__neighbour_found += 1

      if tracks_expected :
        self.__expected_track_events += 1
      else :
        self.__no_expected_track_events += 1

      # if expecting to see a track in the other tracker
      if tracks_expected_other:
        self.__num_expected_other_tracker += 1
        if total_tracks == 0:
          self.__missing_here += 1
      else: #not expecting in other tracker
        self.__no_expected_other += 1
      
      if tracks_expected and total_tracks == 0:
        # expected a track upstream
        self.__expected_missing_track_events += 1

      if has_neighbours and total_tracks == 0:
        self.__neighbour_missed += 1

      if total_tracks == 0 and has_neighbours and tracks_expected :
        self.__neighbour_expected_missed += 1


      for track in kalman_tracks :
        self.__plot_chisq_ndf.Fill(track.chi2()/track.ndf())
        # track points
        tps = track.scifitrackpoints()
        for tp in tps:
            t = tp.tracker()
            s = tp.station()
            p = tp.plane()
            x = tp.pos().x()
            y = tp.pos().y()
            z = tp.pos().z()
            for cl in clusters:
                clt = cl.get_tracker()
                cls = cl.get_station()
                clp = cl.get_plane()
                clused = cl.is_used()
                clx = cl.get_position().x()
                cly = cl.get_position().y()
                clz = cl.get_position().z()


      for key, value in spacepoints_station.iteritems() :
        self.__plot_spacepoints_station_good.Fill( key, value )

      for cluster in clusters :
        if cluster.get_tracker() == self.__tracker :
          s = cluster.get_station()-1
          p = cluster.get_plane()
          self.__plot_cluster_npe.Fill(cluster.get_npe())
          if tracks_expected:
            self.__plot_cluster_npe_expected.Fill(cluster.get_npe())
          else:
            self.__plot_cluster_npe_notexpected.Fill(cluster.get_npe())
          if has_neighbours:
            self.__plot_cluster_npe_neighbours.Fill(cluster.get_npe())
          else:
            self.__plot_cluster_npe_noneighbours.Fill(cluster.get_npe())
          self.__plot_clnpe[s][p].Fill(cluster.get_npe())
          if tracks_expected:
            self.__plot_clnpe_expected[s][p].Fill(cluster.get_npe())
          else:
            self.__plot_clnpe_no_expected[s][p].Fill(cluster.get_npe())


      for digit in digits:
        if digit.get_tracker() == self.__tracker :
          s = digit.get_station()-1
          p = digit.get_plane()
          c = digit.get_channel()
          self.__plot_dignpeVchan[s][p].Fill(c, digit.get_npe())
          self.__plot_dignpe[s][p].Fill(digit.get_npe())
          self.__plot_digadc[s][p].Fill(digit.get_adc())
          self.__plot_digits[s][p].Fill(c)
          if tracks_expected:
              self.__plot_dignpe_expected[s][p].Fill(digit.get_npe())

      for spacepoint in spacepoints :
        if spacepoint.get_tracker() == self.__tracker :
          self.__plot_spacepoint_npe.Fill(spacepoint.get_npe())
          
          spx = spacepoint.get_global_position().X()
          spy = spacepoint.get_global_position().Y()
          stn = spacepoint.get_station()-1
          if not spacepoint.is_used():
            self.__plot_spnpe_notused[stn].Fill( spacepoint.get_npe() )
            self.__plot_spxy_notused[stn].Fill(spx, spy)
          if total_tracks == 1 :
            if spacepoint.is_used() :
              self.__plot_found_spacepoint_npe.Fill( spacepoint.get_npe() )
              self.__plot_spnpe_used[stn].Fill( spacepoint.get_npe() )
              self.__plot_spxy_used[stn].Fill(spx, spy)
            else :
              self.__plot_unused_spacepoint_npe.Fill( spacepoint.get_npe() )
          if tracks_expected:
              self.__plot_expected_spnpe.Fill( spacepoint.get_npe() )
              self.__plot_expected_spxy.Fill(spx, spy)
          else:
              self.__plot_notexpected_spnpe.Fill( spacepoint.get_npe() )
              self.__plot_notexpected_spxy.Fill(spx, spy)
          if has_neighbours:
              self.__plot_neighbours_spnpe.Fill( spacepoint.get_npe() )
              self.__plot_neighbours_spxy.Fill(spx, spy)
          else:
              self.__plot_noneighbours_spnpe.Fill( spacepoint.get_npe() )
              self.__plot_noneighbours_spxy.Fill(spx, spy)

      if num_spacepoints_CUT <= 3:
        for spacepoint in spacepoints :
          stn = spacepoint.get_station()-1
          if spacepoint.get_tracker() == self.__tracker :
            self.__plot_cutevt_spacepoint_npe.Fill(spacepoint.get_npe())
          
            spx = spacepoint.get_global_position().X()
            spy = spacepoint.get_global_position().Y()
            self.__plot_spxy_cutevt[stn].Fill(spx, spy)
            self.__plot_spnpe_cutevt[stn].Fill( spacepoint.get_npe() )
            if tracks_expected:
              self.__plot_cutevt_expected_spnpe.Fill( spacepoint.get_npe() )
              self.__plot_cutevt_expected_spxy.Fill(spx, spy)
            else:
              self.__plot_cutevt_notexpected_spnpe.Fill( spacepoint.get_npe() )
              self.__plot_cutevt_notexpected_spxy.Fill(spx, spy)
            if has_neighbours:
              self.__plot_cutevt_neighbours_spnpe.Fill( spacepoint.get_npe() )
              self.__plot_cutevt_neighbours_spxy.Fill(spx, spy)
            else:
              self.__plot_cutevt_noneighbours_spnpe.Fill( spacepoint.get_npe() )
              self.__plot_cutevt_noneighbours_spxy.Fill(spx, spy)

      if num_spacepoints_CUT > 3 :
        self.__good_event_counter += 1
        if total_tracks < 1 :
          self.__old_select_missed += 1
      else :
        if total_tracks > 0 :
          self.__new_select_missed += 1

      if num_helicals > 0 :
        self.__helical_events += 1
      if num_straights > 0 :
        self.__straight_events += 1

      self.__plot_spacepoints.Fill(num_spacepoints)
      self.__plot_spacepoints_CUT.Fill(num_spacepoints_CUT)

    else :
      self.__bad_event_counter += 1

      if num_spacepoints_CUT > 2 :
        if total_tracks < 1 :
          self.__old_missed += 1
        else :
          self.__old_found += 1
      else :
        if total_tracks > 0 :
          self.__new_missed += 1
        else :
          self.__new_found += 1


  ###############################################
  def conclude(self) :
    print
    print "Tracker Number : {0}".format(self.__tracker)
    print
    print "Analysed {0} Events".format(self.__event_counter)
    print 
    print "Spacepoint NPE Cut {0}:  ".format(NPE_CUT)
    print 
    print "Mean Numbers :"
    print "  Digits      :", float(self.__digit_counter)/float(self.__event_counter)
    print "  Clusters    :", float(self.__cluster_counter)/float(self.__event_counter)
    print "  Spacepoints :", float(self.__spacepoint_counter)/float(self.__event_counter)
    print "  Helicals    :", float(self.__helical_counter)/float(self.__event_counter)
    print "  Straights   :", float(self.__straight_counter)/float(self.__event_counter)
    print
    print "  Spacepoints > {0:3d}        : {1:2.4f}".format( int(NPE_CUT), float(self.__spacepoint_CUT_counter)/float(self.__event_counter) )
    print "  Spacepoints > {0:3d} (Good) : {1:2.4f}".format( int(NPE_CUT), float(self.__spacepoint_CUT_counter)/float(self.__good_event_counter) )
    print
    print "Total Numbers :"
    print "  Digits      :", self.__digit_counter
    print "  Clusters    :", self.__cluster_counter
    print "  Spacepoints :", self.__spacepoint_counter
    print "  Helicals    :", self.__helical_counter
    print "  Straights   :", self.__straight_counter
    print
    print "  Events with 3xSpacepoints > {0:3d} : {1:3d}".format( int(NPE_CUT), self.__good_event_counter )
    print
    print "  Events with at least 1 helical  :", self.__helical_events
    print "  Events with at least 1 straight :", self.__straight_events
    print
    print "New Selection Dropped :", self.__new_select_missed
    print "Old Selection Dropped :", self.__old_select_missed
    print
    print " New Selection Found/Dropped: ", self.__new_found, self.__new_missed
    print " Old Selection Found/Dropped: ", self.__old_found, self.__old_missed
    print
    print "Total events before cuts:         ", self.__total_event_counter
    print "Events that passed cuts:          ", self.__event_counter
    print "Events that didn't pass the cuts :", self.__bad_event_counter
    print
    print "Events Expecting a track:         ", self.__expected_track_events
    print "Events Not Expecting a track:     ", self.__no_expected_track_events
    print "Events Missing an expected track: ", self.__expected_missing_track_events
    print
    print "Events expecting in other tracker:         ", self.__num_expected_other_tracker
    print "Track in other tracker, but missing here:  ", self.__missing_here
    print "Events Not expecting in other tracker:     ", self.__no_expected_other
    print
    print "Events with a neighbour                 : ", self.__neighbour_found
    print "Events with a neighbour, missing a track: ", self.__neighbour_missed, " = ", 100.0*(self.__neighbour_found - self.__neighbour_missed) / self.__neighbour_found, "%"
    print "Events with a neighbour, and expected, missing a track: ", self.__neighbour_expected_missed, " = ", 100.0*(self.__neighbour_found - self.__neighbour_expected_missed) / self.__neighbour_found, "%"
    print
    print


  ###############################################
  def get_plots(self) :
    plot_dict = {}

    for s in range(5):
        name = "sp_npe_used_"+str(s+1)
        plot_dict[name] = self.__plot_spnpe_used[s]
        name = "sp_npe_notused_"+str(s+1)
        plot_dict[name] = self.__plot_spnpe_notused[s]
        name = "sp_xy_used_"+str(s+1)
        plot_dict[name] = self.__plot_spxy_used[s]
        name = "sp_xy_notused_"+str(s+1)
        plot_dict[name] = self.__plot_spxy_notused[s]
        name = "sp_npe_cutevt"+str(s+1)
        plot_dict[name] = self.__plot_spnpe_cutevt[s]
        name = "sp_xy_cutevt"+str(s+1)
        plot_dict[name] = self.__plot_spxy_cutevt[s]
        for p in range(3):
            name = "cluster_npe_s"+str(s+1)+"_p"+str(p)
            plot_dict[name] = self.__plot_clnpe[s][p]
            name = "cluster_npe_expected_s"+str(s+1)+"_p"+str(p)
            plot_dict[name] = self.__plot_clnpe_expected[s][p]
            name = "cluster_npe_not_expected_s"+str(s+1)+"_p"+str(p)
            plot_dict[name] = self.__plot_clnpe_no_expected[s][p]
            name = "digit_npe_s"+str(s+1)+"_p"+str(p)
            plot_dict[name] = self.__plot_dignpe[s][p]
            name = "digit_npe_expected_s"+str(s+1)+"_p"+str(p)
            plot_dict[name] = self.__plot_dignpe_expected[s][p]
            name = "digit_npeVchan_s"+str(s+1)+"_p"+str(p)
            plot_dict[name] = self.__plot_dignpeVchan[s][p]
            name = "digit_adc_s"+str(s+1)+"_p"+str(p)
            plot_dict[name] = self.__plot_digadc[s][p]
            name = "digits_"+str(s+1)+"_p"+str(p)
            plot_dict[name] = self.__plot_digits[s][p]
    plot_dict['cluster_npe'] = self.__plot_cluster_npe
    plot_dict['cluster_npe_expected'] = self.__plot_cluster_npe_expected
    plot_dict['cluster_npe_notexpected'] = self.__plot_cluster_npe_notexpected
    plot_dict['cluster_npe_neighbours'] = self.__plot_cluster_npe_neighbours
    plot_dict['cluster_npe_noneighbours'] = self.__plot_cluster_npe_noneighbours
    plot_dict['spacepoint_npe'] = self.__plot_spacepoint_npe
    plot_dict['cutevt_spacepoint_npe'] = self.__plot_cutevt_spacepoint_npe
    plot_dict['spacepoints'] = self.__plot_spacepoints
    plot_dict['spacepoints_CUT'] = self.__plot_spacepoints_CUT

    plot_dict['found_spacepoint_npe'] = self.__plot_found_spacepoint_npe
    plot_dict['unused_spacepoint_npe'] = self.__plot_unused_spacepoint_npe

    plot_dict['expected_spnpe'] = self.__plot_expected_spnpe
    plot_dict['notexpected_spnpe'] = self.__plot_notexpected_spnpe
    plot_dict['expected_spxy'] = self.__plot_expected_spxy
    plot_dict['notexpected_spxy'] = self.__plot_notexpected_spxy

    plot_dict['cutevt_expected_spnpe'] = self.__plot_cutevt_expected_spnpe
    plot_dict['cutevt_notexpected_spnpe'] = self.__plot_cutevt_notexpected_spnpe
    plot_dict['cutevt_expected_spxy'] = self.__plot_cutevt_expected_spxy
    plot_dict['cutevt_notexpected_spxy'] = self.__plot_cutevt_notexpected_spxy

    plot_dict['neighbours_spnpe'] = self.__plot_neighbours_spnpe
    plot_dict['noneighbours_spnpe'] = self.__plot_noneighbours_spnpe
    plot_dict['neighbours_spxy'] = self.__plot_neighbours_spxy
    plot_dict['noneighbours_spxy'] = self.__plot_noneighbours_spxy

    plot_dict['cutevt_neighbours_spnpe'] = self.__plot_cutevt_neighbours_spnpe
    plot_dict['cutevt_noneighbours_spnpe'] = self.__plot_cutevt_noneighbours_spnpe
    plot_dict['cutevt_neighbours_spxy'] = self.__plot_cutevt_neighbours_spxy
    plot_dict['cutevt_noneighbours_spxy'] = self.__plot_cutevt_noneighbours_spxy

    plot_dict['new_spacepoints_CUT'] = self.__plot_spacepoints_CUT_new
    plot_dict['old_spacepoints_CUT'] = self.__plot_spacepoints_CUT_old

    plot_dict['spacepoints_station_good'] = self.__plot_spacepoints_station_good
    plot_dict['spacepoints_station_new'] = self.__plot_spacepoints_station_new
    plot_dict['spacepoints_station_old'] = self.__plot_spacepoints_station_old

    return plot_dict



if __name__ == "__main__" : 
  ROOT.gROOT.SetBatch( True )
  ROOT.gErrorIgnoreLevel = ROOT.kError

  parser = argparse.ArgumentParser( description="" )

  parser.add_argument( 'maus_root_files', nargs='+', help='List of MAUS '+\
                  'output root files containing reconstructed straight tracks')

  try :
    namespace = parser.parse_args()

    upstream_analyser = Tracker_Analyser(0)
    downstream_analyser = Tracker_Analyser(1)

  except BaseException as ex:
    raise
  else :


    file_reader = event_loader.maus_reader(namespace.maus_root_files)
    file_reader.set_max_num_events(50000)

    try :
      while file_reader.next_selected_event() :

        try :
          scifi_event = file_reader.get_event( 'scifi' )
          tof_event = file_reader.get_event( 'tof' )

          if passes_cuts( tof_event ) :
            upstream_analyser.analyse_event( scifi_event, True )
            downstream_analyser.analyse_event( scifi_event, True )
          else :
            upstream_analyser.analyse_event( scifi_event, False )
            downstream_analyser.analyse_event( scifi_event, False )

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


    print"\n- Analysing Data...\n"

    upstream_analyser.conclude()
    downstream_analyser.conclude()

    plots = { 'upstream' : upstream_analyser.get_plots(), 'downstream' : downstream_analyser.get_plots(), 'tof01' : TOF01_PLOT, 'tof12' : TOF12_PLOT }

    analysis.tools.save_plots(plots, "./pr_analysis.root")
    
