#include "JEventProcessor_EVIO.h"
#include <JANA/JLogger.h>

/**
 * @brief Constructor for JEventProcessor_EVIO
 * 
 * Initialize the processor with the appropriate type name, prefix, and callback style.
 */
JEventProcessor_EVIO::JEventProcessor_EVIO() {
    SetTypeName(NAME_OF_THIS);                    // Provide JANA with this class's name
    SetPrefix("jeventprocessor_evio");            // Set unique prefix for parameters
    SetCallbackStyle(CallbackStyle::ExpertMode);  // Use expert mode for full control

    // All of these are optional because not all events will have these hits
    m_caen1190_hits_in.SetOptional(true);
    m_fadc_scaler_hits_in.SetOptional(true);
    m_pulse_hits_in.SetOptional(true);
    m_waveform_hits_in.SetOptional(true);
    m_ti_scaler_hits_in.SetOptional(true);
    m_heldec_data_in.SetOptional(true);
    m_mpd_hits_in.SetOptional(true);
    m_vftdc_hits_in.SetOptional(true);
    m_hallb_pulse_integral_hits_in.SetOptional(true);
    m_hallb_pulse_time_hits_in.SetOptional(true);
    m_hallb_pulse_peak_hits_in.SetOptional(true);
    m_faV3compton_hits_in.SetOptional(true);
    m_faV3comptonAccumulator_hits_in.SetOptional(true);
}

/**
 * @brief Initialize the processor
 * 
 * Called once at the start of processing. Open the output files and set up
 * any necessary resources for event processing.
 */
void JEventProcessor_EVIO::Init() {
    LOG << "JEventProcessor_EVIO::Init" << LOG_END;
    
    // Open the ROOT output file
    m_root_output_file = new TFile(m_root_output_filename().c_str(), "RECREATE");
    if (m_root_output_file == nullptr || m_root_output_file->IsZombie()) {
        throw JException("Failed to open ROOT output file: " + m_root_output_filename());  
    }

    // Initialize the waveform tree row data structure
    m_waveform_tree_row = WaveformTreeRow();

    // Create ROOT tree for waveform data
    m_waveform_tree = new TTree("waveform_tree", "FADC250 Waveform Data (slot, channel, waveform)");
    m_waveform_tree->Branch("slot", &ev_slot);
    m_waveform_tree->Branch("chan", &ev_chan);
    m_waveform_tree->Branch("waveform", &ev_waveform);
    m_waveform_tree->Branch("rocid", &ev_rocid);

    // Create ROOT tree for pulse data
    m_pulse_tree = new TTree("pulse_tree","FADC250 pulse data(slow, channel, integral, time)");
    m_pulse_tree->Branch("integral_sum", &ev_integral_sum);
    m_pulse_tree->Branch("pedestal_sum", &pedestal_sum);
    m_pulse_tree->Branch("coarse_time",&ev_coarse_time);
    m_pulse_tree->Branch("fine_time",&ev_fine_time);
    m_pulse_tree->Branch("pulse_peak",&ev_pulse_peak);
    m_pulse_tree->Branch("pedestal_quality",&pedestal_quality);
    m_pulse_tree->Branch("nhits",&number_hit);
    m_pulse_tree->Branch("chan",&ev_pulse_chan);
    m_pulse_tree->Branch("slot",&ev_pulse_slot);
    m_pulse_tree->Branch("rocid", &ev_pulse_rocid);

    // Create ROOT tree for CAEN1190 data
    m_caen1190_tree = new TTree("caen1190_tree","CAEN1190 data (slot, channel, time)");
    m_caen1190_tree->Branch("rocid", &ev_caen_rocid);
    m_caen1190_tree->Branch("slot", &ev_caen_slot);
    m_caen1190_tree->Branch("chan", &ev_caen_chan);
    m_caen1190_tree->Branch("measurement", &ev_caen_measurement);
    m_caen1190_tree->Branch("opt", &ev_caen_opt);
    m_caen1190_tree->Branch("flags", &ev_caen_flags);
    m_caen1190_tree->Branch("trig_time", &ev_caen_trig_time);
    m_caen1190_tree->Branch("hdr_chip_id", &ev_caen_hdr_chip_id);
    m_caen1190_tree->Branch("hdr_event_id", &ev_caen_hdr_event_id);
    m_caen1190_tree->Branch("hdr_bunch_id", &ev_caen_hdr_bunch_id);
    m_caen1190_tree->Branch("trl_status", &ev_caen_trl_status);

    // Create the Helicity Decoder Tree
    m_tree = new TTree("m_tree", "Physics Event Tree");
    m_tree->Branch(
        "heldec",
        &heldec,
        "helicity_seed/i:"
        "n_tstable_fall/i:"
        "n_tstable_rise/i:"
        "n_pattsync/i:"
        "n_pairsync/i:"
        "time_tstable_start/i:"
        "time_tstable_end/i:"
        "last_tstable_duration/i:"
        "last_tsettle_duration/i:"
        "trig_tstable/i:"
        "trig_pattsync/i:"
        "trig_pairsync/i:"
        "trig_helicity/i:"
        "trig_pat0_helicity/i:"
        "trig_polarity/i:"
        "trig_pat_count/i:"
        "last32wins_pattsync/i:"
        "last32wins_pairsync/i:"
        "last32wins_helicity/i:"
        "last32wins_pattsync_hel/i"
    );

    // faV3 compton tree
    compton_tree = new TTree("compton_tree", "compton tree");
    compton_tree->Branch("event_number", &comp_trigger_num);
    compton_tree->Branch("trig_timestamp1", &comp_timestamp1);
    compton_tree->Branch("trig_timestamp2", &comp_timestamp2);
    compton_tree->Branch("rocid", &comp_rocid);
    compton_tree->Branch("slot", &comp_slot);
    compton_tree->Branch("moduleid", &comp_module_id);

    compton_tree->Branch("acc_samp_overflow", &acc_samp_overflow);
    compton_tree->Branch("acc_samp_underflow", &acc_samp_underflow);
    compton_tree->Branch("acc_type", &acc_type);
    compton_tree->Branch("acc_chan", &acc_chan);
    compton_tree->Branch("acc_overflow_timestamp1", &acc_overflow_timestamp1);
    compton_tree->Branch("acc_overflow_timestamp2", &acc_overflow_timestamp2);
    compton_tree->Branch("acc_underflow_timestamp1", &acc_underflow_timestamp1);
    compton_tree->Branch("acc_underflow_timestamp2", &acc_underflow_timestamp2);
    compton_tree->Branch("acc_sum_nsample1", &acc_sum_nsample1);
    compton_tree->Branch("acc_sum_nsample2", &acc_sum_nsample2);
    compton_tree->Branch("acc_sum1", &acc_sum1);
    compton_tree->Branch("acc_sum2", &acc_sum2);
    compton_tree->Branch("acc_sum", &acc_sum);
    compton_tree->Branch("acc_sum_nsample", &acc_sum_nsample);
    compton_tree->Branch("acc_np_nsboverlapped", &acc_np_nsboverlapped);
    compton_tree->Branch("acc_np_nonsa", &acc_np_nonsa);
    compton_tree->Branch("acc_np_miss", &acc_np_miss);

    // Optionally: Text output file for human-readable hit summaries
    m_txt_output_file.open(m_txt_output_filename().c_str());

    // Create histogram for pulse integral distribution
    m_pulse_integral_hist = new TH1I("h_integral", "Pulse Integral Distribution;Integral Sum;Counts", 100, 0, 1);
    m_pulse_integral_hist->SetCanExtend(TH1::kAllAxes);  // Allow ROOT to automatically extend bins
}

/**
 * @brief Process a single event sequentially
 * 
 * Processes FADC250 detector data for a single event. Fills ROOT tree with
 * waveform data and histogram with pulse integral values. This method is 
 * called for each event in the processing pipeline.
 * 
 * @param event Reference to the JANA2 event to process
 */
void JEventProcessor_EVIO::ProcessSequential(const JEvent &event) {
    
    // Clear previous event data
    ev_slot.clear();
    ev_chan.clear();
    ev_waveform.clear();
    ev_rocid.clear();
    ev_coarse_time.clear();
    ev_pulse_chan.clear();
    ev_pulse_slot.clear();
    ev_fine_time.clear();
    ev_integral_sum.clear();
    ev_pulse_peak.clear();
    ev_pulse_rocid.clear();
    pedestal_sum= 0;
    pedestal_quality = 0;
    number_hit =0;

    // Clear previous event data - CAEN1190
    ev_caen_rocid.clear();
    ev_caen_slot.clear();
    ev_caen_chan.clear();
    ev_caen_measurement.clear();
    ev_caen_opt.clear();
    ev_caen_flags.clear();
    ev_caen_trig_time.clear();
    ev_caen_hdr_chip_id.clear();
    ev_caen_hdr_event_id.clear();
    ev_caen_hdr_bunch_id.clear();
    ev_caen_trl_status.clear();

    // CAEN1190 TDC hits
    for (const auto& caen_hit : m_caen1190_hits_in()) {
        ev_caen_rocid.push_back(caen_hit->rocid);
        ev_caen_slot.push_back(caen_hit->slot);
        ev_caen_chan.push_back(caen_hit->chan);
        ev_caen_measurement.push_back(caen_hit->measurement);
        ev_caen_opt.push_back(caen_hit->opt);
        ev_caen_flags.push_back(caen_hit->flags);
        ev_caen_trig_time.push_back(caen_hit->trig_time);
        ev_caen_hdr_chip_id.push_back(caen_hit->hdr_chip_id);
        ev_caen_hdr_event_id.push_back(caen_hit->hdr_event_id);
        ev_caen_hdr_bunch_id.push_back(caen_hit->hdr_bunch_id);
        ev_caen_trl_status.push_back(caen_hit->glb_trl_status);

    }
    m_caen1190_tree->Fill();

    // FADC250 waveform hits
    for (const auto& waveform_hit : m_waveform_hits_in()) {
        // Fill ROOT tree with waveform data
        m_waveform_tree_row.slot = waveform_hit->slot;
        m_waveform_tree_row.chan = waveform_hit->chan;
        m_waveform_tree_row.rocid = waveform_hit->rocid;
        m_waveform_tree_row.waveform = waveform_hit->waveform;

	size_t waveform_sample_number = m_waveform_tree_row.waveform.size();

	ev_slot.insert(ev_slot.end(), waveform_sample_number, m_waveform_tree_row.slot);
	ev_chan.insert(ev_chan.end(), waveform_sample_number, m_waveform_tree_row.chan);
        ev_rocid.insert(ev_rocid.end(), waveform_sample_number, m_waveform_tree_row.rocid);
	ev_waveform.insert(ev_waveform.end(),  m_waveform_tree_row.waveform.begin(), m_waveform_tree_row.waveform.end());
    }

    // FADC250 pulse hits
 
    int nn=0;
    for (const auto& pulse_hit : m_pulse_hits_in()){
        integral_sum = pulse_hit->integral_sum;
        pedestal_sum = pulse_hit->pedestal_sum;
        coarse_time = pulse_hit->coarse_time;
        fine_time = pulse_hit->fine_time;
        pulse_peak = pulse_hit->pulse_peak;
        pedestal_quality = pulse_hit->pedestal_quality;
        if(integral_sum!=0){
            nn++;
            ev_integral_sum.push_back(integral_sum);
            ev_coarse_time.push_back(coarse_time);
            ev_fine_time.push_back(fine_time);
            ev_pulse_peak.push_back(pulse_peak);
            ev_pulse_slot.push_back(pulse_hit->slot);
            ev_pulse_chan.push_back(pulse_hit->chan);
            ev_pulse_rocid.push_back(pulse_hit->rocid);
        }
       
        
    }

    number_hit = nn;
    m_waveform_tree->Fill();
    if(nn>0){
        m_pulse_tree->Fill();
    }


    // FADC250 pulse hits
    for (const auto& pulse_hit : m_pulse_hits_in()) {
        // Fill histogram with pulse integral values
        m_pulse_integral_hist->Fill(pulse_hit->integral_sum);
    }

    heldec = {};
    // Helicity decoder data
    for(const auto& heldec_hit : m_heldec_data_in()){
	heldec.helicity_seed          = heldec_hit->helicity_seed;
        heldec.n_tstable_fall         = heldec_hit->n_tstable_fall;
        heldec.n_tstable_rise         = heldec_hit->n_tstable_rise;
        heldec.n_pattsync             = heldec_hit->n_pattsync;
        heldec.n_pairsync             = heldec_hit->n_pairsync;
        heldec.time_tstable_start     = heldec_hit->time_tstable_start;
        heldec.time_tstable_end       = heldec_hit->time_tstable_end;
        heldec.last_tstable_duration  = heldec_hit->last_tstable_duration;
        heldec.last_tsettle_duration  = heldec_hit->last_tsettle_duration;
        heldec.trig_tstable           = heldec_hit->trig_tstable;
        heldec.trig_pattsync          = heldec_hit->trig_pattsync;
        heldec.trig_pairsync          = heldec_hit->trig_pairsync;
        heldec.trig_helicity          = heldec_hit->trig_helicity;
        heldec.trig_pat0_helicity     = heldec_hit->trig_pat0_helicity;
        heldec.trig_polarity          = heldec_hit->trig_polarity;
        heldec.trig_pat_count         = heldec_hit->trig_pat_count;
        heldec.last32wins_pattsync    = heldec_hit->last32wins_pattsync;
        heldec.last32wins_pairsync    = heldec_hit->last32wins_pairsync;
        heldec.last32wins_helicity    = heldec_hit->last32wins_helicity;
        heldec.last32wins_pattsync_hel= heldec_hit->last32wins_pattsync_hel;
        m_tree->Fill();
    }

    for(const auto& compton_evt : m_faV3compton_hits_in()){
        comp_trigger_num = compton_evt->trigger_num;
	comp_timestamp1 = compton_evt->timestamp1;
	comp_timestamp2 = compton_evt->timestamp2;
        comp_rocid = compton_evt->rocid;
	comp_slot = compton_evt->slot;
	comp_module_id = compton_evt->module_id;
    }

    acc_samp_overflow.clear();
    acc_samp_underflow.clear();
    acc_type.clear();
    acc_chan.clear();
    acc_overflow_timestamp1.clear();
    acc_overflow_timestamp2.clear();
    acc_underflow_timestamp1.clear();
    acc_underflow_timestamp2.clear();
    acc_sum_nsample1.clear();
    acc_sum_nsample2.clear();
    acc_sum1.clear();
    acc_sum2.clear();
    acc_sum.clear();
    acc_sum_nsample.clear();
    acc_np_nsboverlapped.clear();
    acc_np_nonsa.clear();
    acc_np_miss.clear();
    for(const auto& compton_accum : m_faV3comptonAccumulator_hits_in()){
        acc_samp_overflow.push_back(compton_accum->acc_samp_overflow);
        acc_samp_underflow.push_back(compton_accum->acc_samp_underflow);
	acc_type.push_back(compton_accum->acc_type);
	acc_chan.push_back(compton_accum->acc_chan);
	acc_overflow_timestamp1.push_back(compton_accum->acc_overflow_timestamp1);
	acc_overflow_timestamp2.push_back(compton_accum->acc_overflow_timestamp2);
	acc_underflow_timestamp1.push_back(compton_accum->acc_underflow_timestamp1);
	acc_underflow_timestamp2.push_back(compton_accum->acc_underflow_timestamp2);
	acc_sum_nsample1.push_back(compton_accum->acc_sum_nsample1);
	acc_sum_nsample2.push_back(compton_accum->acc_sum_nsample2);
	acc_sum1.push_back(compton_accum->acc_sum1);
	acc_sum2.push_back(compton_accum->acc_sum2);
    acc_sum.push_back(compton_accum->acc_sum);
    acc_sum_nsample.push_back(compton_accum->acc_sum_nsample);
	acc_np_nsboverlapped.push_back(compton_accum->acc_np_nsboverlapped);
	acc_np_nonsa.push_back(compton_accum->acc_np_nonsa);
	acc_np_miss.push_back(compton_accum->acc_np_miss);
    }

    compton_tree->Fill();

    // ------------------------------------------------------------------
    // Optional text dump of hits for this event (waveforms, pulses, scalers)
    // ------------------------------------------------------------------
    if (m_txt_output_file.is_open()) {
        const auto& caen1190_hits            = m_caen1190_hits_in();
        const auto& waveform_hits            = m_waveform_hits_in();
        const auto& pulse_hits               = m_pulse_hits_in();
        const auto& fadc_scaler_hits         = m_fadc_scaler_hits_in();
        const auto& ti_scaler_hits           = m_ti_scaler_hits_in();
        const auto& mpd_hits                 = m_mpd_hits_in();
        const auto& vftdc_hits               = m_vftdc_hits_in();
        const auto& hallb_pulse_integral_hits = m_hallb_pulse_integral_hits_in();
        const auto& hallb_pulse_time_hits    = m_hallb_pulse_time_hits_in();
        const auto& hallb_pulse_peak_hits    = m_hallb_pulse_peak_hits_in();
        const auto& faV3compton_hits    = m_faV3compton_hits_in();

        bool have_caen1190_hits          = !caen1190_hits.empty();
        bool have_waveforms              = !waveform_hits.empty();
        bool have_pulses                 = !pulse_hits.empty();
        bool have_fadc_scalers           = !fadc_scaler_hits.empty();
        bool have_ti_scalers             = !ti_scaler_hits.empty();
        bool have_mpd_hits               = !mpd_hits.empty();
        bool have_vftdc_hits             = !vftdc_hits.empty();
        bool have_hallb_pulse_integrals  = !hallb_pulse_integral_hits.empty();
        bool have_hallb_pulse_times      = !hallb_pulse_time_hits.empty();
        bool have_hallb_pulse_peaks      = !hallb_pulse_peak_hits.empty();
        // Only write anything if we have at least one type of hit
        if (have_caen1190_hits ||  have_waveforms || have_pulses || have_fadc_scalers || have_ti_scalers || have_mpd_hits || have_vftdc_hits
            || have_hallb_pulse_integrals || have_hallb_pulse_times || have_hallb_pulse_peaks) {
            auto event_number = event.GetEventNumber();

            m_txt_output_file << "Event " << event_number << "\n";

            // CAEN1190 summary
            if (have_caen1190_hits) {
                m_txt_output_file << "  CAEN1190 hits: " << caen1190_hits.size() << "\n";
                for (const auto& hit : caen1190_hits) {
                    m_txt_output_file
                        << "    CAEN1190 rocid=" << hit->rocid
                        << " slot=" << hit->slot
                        << " chan=" << hit->chan
                        << " measurement=" << hit->measurement
                        << " opt=" << hit->opt
                        << " flags=" << hit->flags
                        << " trig_time=" << hit->trig_time
                        << " hdr_chip_id=" << hit->hdr_chip_id
                        << " hdr_event_id=" << hit->hdr_event_id
                        << " hdr_bunch_id=" << hit->hdr_bunch_id
                        << " trl_status=" << hit->glb_trl_status
                        << "\n";
                }
            } else {
                m_txt_output_file << "  No CAEN1190hit in this event\n";
            }

            // Waveform summary
            if (have_waveforms) {
                m_txt_output_file << "  Waveform hits: " << waveform_hits.size() << "\n";
                for (const auto& hit : waveform_hits) {
                    m_txt_output_file
                        << "    WF rocid=" << hit->rocid
                        << " slot=" << hit->slot
                        << " chan=" << hit->chan
                        << " nsamples=" << hit->waveform.size()
                        << "\n";
                }
            } else {
                m_txt_output_file << "  No FADC250 waveform hits in this event\n";
            }

            // Pulse summary
            if (have_pulses) {
                m_txt_output_file << "  Pulse hits: " << pulse_hits.size() << "\n";
                for (const auto& hit : pulse_hits) {
                    m_txt_output_file
                        << "    PULSE rocid=" << hit->rocid
                        << " slot=" << hit->slot
                        << " chan=" << hit->chan
                        << " integral_sum=" << hit->integral_sum
                        << "\n";
                }
            } else {
                m_txt_output_file << "  No FADC250 pulse hits in this event\n";
            }

            // FADC scaler summary
            if (have_fadc_scalers) {
                m_txt_output_file << "  FADC scaler hits: " << fadc_scaler_hits.size() << "\n";
                for (const auto& hit : fadc_scaler_hits) {
                    m_txt_output_file
                        << "    SCALER rocid=" << hit->rocid
                        << " slot=" << hit->slot
                        << " ncounts=" << hit->ncounts
                        << " counts=";
                    for (uint32_t i = 0; i < hit->ncounts && i < 16u; ++i) {
                        m_txt_output_file << hit->counts[i];
                        if (i + 1 < hit->ncounts && i + 1 < 16u) {
                            m_txt_output_file << ",";
                        }
                    }
                    m_txt_output_file << "\n";
                }
            } else {
                m_txt_output_file << "  No FADCScalerHit objects in this event\n";
            }

            // TI scaler summary
            if (have_ti_scalers) {
                m_txt_output_file << "  TI scaler hits: " << ti_scaler_hits.size() << "\n";
                for (const auto& hit : ti_scaler_hits) {
                    m_txt_output_file
                        << "    TISCALER rocid=" << hit->rocid
                        << " slot=" << hit->slot
                        << " nwords=" << hit->nscalerwords
                        << " live_time=" << hit->live_time
                        << " busy_time=" << hit->busy_time
                        << " ts_inputs_before_busy=" << hit->ts_inputs_before_busy
                        << "\n";
                }
            } else {
                m_txt_output_file << "  No TIScalerHit objects in this event\n";
            }

            // MPD hit summary
            if (have_mpd_hits) {
                m_txt_output_file << "  MPD hits: " << mpd_hits.size() << "\n";
                for (const auto& hit : mpd_hits) {
                    m_txt_output_file
                        << "    MPD rocid=" << hit->rocid
                        << " slot=" << hit->slot
                        << " trigger_num=" << hit->trigger_num
                        << " trigger_time=" << hit->trigger_time
                        << " mpd_id=" << (int)hit->mpd_id
                        << " fiber_id=" << (int)hit->fiber_id
                        << " apv_channel=" << (int)hit->apv_channel
                        << " apv_id=" << (int)hit->apv_id
                        << " apv_samples=[";
                    for (int i = 0; i < 6; ++i) {
                        m_txt_output_file << hit->apv_samples[i];
                        if (i < 5) m_txt_output_file << ",";
                    }
                    m_txt_output_file << "]\n";
                }
            } else {
                m_txt_output_file << "  No MPDHit objects in this event\n";
            }

            // VFTDC hit summary
            if (have_vftdc_hits) {
                m_txt_output_file << "  VFTDC hits: " << vftdc_hits.size() << "\n";
                for (const auto& hit : vftdc_hits) {
                    m_txt_output_file
                        << "    VFTDC rocid=" << hit->rocid
                        << " slot=" << hit->slot
                        << " board_id=" << hit->board_id
                        << " timestamp=" << hit->timestamp
                        << " group_num=" << hit->group_num
                        << " channel_num=" << hit->channel_num
                        << " edge_type=" << hit->edge_type
                        << " coarse_time=" << hit->coarse_time
                        << " fine_time=" << hit->fine_time
                        << " two_ns=" << hit->two_ns
                        << "\n";
                }
            } else {
                m_txt_output_file << "  No VFTDCHit objects in this event\n";
            }

            // HallB pulse integral summary
            if (have_hallb_pulse_integrals) {
                m_txt_output_file << "  HallB pulse integral hits: " << hallb_pulse_integral_hits.size() << "\n";
                for (const auto& hit : hallb_pulse_integral_hits) {
                    m_txt_output_file
                        << "    HALLB_INTEGRAL slot=" << hit->slot
                        << " chan=" << hit->chan
                        << " pulse_number=" << hit->pulse_number
                        << " pulse_integral=" << hit->pulse_integral
                        << "\n";
                }
            } else {
                m_txt_output_file << "  No HallB pulse integral hits in this event\n";
            }

            // HallB pulse time summary
            if (have_hallb_pulse_times) {
                m_txt_output_file << "  HallB pulse time hits: " << hallb_pulse_time_hits.size() << "\n";
                for (const auto& hit : hallb_pulse_time_hits) {
                    m_txt_output_file
                        << "    HALLB_TIME slot=" << hit->slot
                        << " chan=" << hit->chan
                        << " pulse_number=" << hit->pulse_number
                        << " measurement_quality_factor=" << hit->measurement_quality_factor
                        << " coarse_pulse_time=" << hit->coarse_pulse_time
                        << " fine_pulse_time=" << hit->fine_pulse_time
                        << "\n";
                }
            } else {
                m_txt_output_file << "  No HallB pulse time hits in this event\n";
            }

            // HallB pulse peak summary
            if (have_hallb_pulse_peaks) {
                m_txt_output_file << "  HallB pulse peak hits: " << hallb_pulse_peak_hits.size() << "\n";
                for (const auto& hit : hallb_pulse_peak_hits) {
                    m_txt_output_file
                        << "    HALLB_PEAK slot=" << hit->slot
                        << " chan=" << hit->chan
                        << " pulse_number=" << hit->pulse_number
                        << " Vmin=" << hit->Vmin
                        << " Vpeak=" << hit->Vpeak
                        << "\n";
                }
            } else {
                m_txt_output_file << "  No HallB pulse peak hits in this event\n";
            }

            m_txt_output_file << "\n";
        }
    }
}

/**
 * @brief Finish processing and cleanup
 * 
 * Called once at the end of processing. Close the output file and perform
 * any necessary cleanup operations.
 */
void JEventProcessor_EVIO::Finish() {
    LOG << "JEventProcessor_EVIO::Finish" << LOG_END;

    // Write ROOT objects and close ROOT file
    if (m_root_output_file) {
        m_waveform_tree->Write();        // Save waveform tree to file
        m_pulse_integral_hist->Write();  // Save integral histogram to file
	    m_tree->Write();
        m_pulse_tree->Write();           // Save pulse tree to file
        m_caen1190_tree->Write();        // Save caen1190 tree to file
        compton_tree->Write();        // Save compton tree to file
        m_root_output_file->Close();     // Close ROOT file
        delete m_root_output_file;       // Free memory
        m_root_output_file = nullptr;
    }

    // Close text output file if open
    if (m_txt_output_file.is_open()) {
        m_txt_output_file.close();
    }
}

