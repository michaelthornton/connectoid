function export_connectoid_results(results, outputDirectory)
% Write the standard CSV set for a spontaneous scan.
arguments
    results struct
    outputDirectory (1,1) string
end
write = @(T,name) writetable(T, fullfile(outputDirectory, name));

write(results.connectionTable, 'connections.csv');
write(results.connectionSummaryTable, 'monosynaptic_connection_counts.csv');
write(results.pairSummaryTable, 'all_directed_pair_summary.csv');
write(results.pairCountTable, 'all_directed_pair_counts.csv');
write(removevars(results.electrodeData,'spikeTimesSec'), 'electrodes.csv');
write(unitExportTable(results.unitData), 'units.csv');
write(roiExportTable(results.organoidRois), 'organoid_rois.csv');

write(results.activity.summaryTable, 'organoid_activity_metrics.csv');
write(results.activity.unitSummaryTable, 'unit_firing_rates.csv');
write(results.activity.peakTable, 'detected_mfr_peaks.csv');

S = results.synchrony;
write(S.summaryTable, 'organoid_synchrony_metrics.csv');
write(S.directionalSummaryTable, 'organoid_synchrony_directional_metrics.csv');
write(S.postTriggerHistogramTable, 'organoid_synchrony_post_trigger_histograms.csv');
write(S.signedHistogramTable, 'organoid_synchrony_signed_histogram.csv');
write(S.postTriggerPairTable, 'organoid_synchrony_peak_pairs.csv');
end

function exported = unitExportTable(unitData)
memberElectrodes = strings(height(unitData),1);
for iUnit = 1:height(unitData)
    memberElectrodes(iUnit) = mat2str(unitData.memberElectrodes{iUnit});
end
exported = table(unitData.unitId, memberElectrodes, unitData.centroidX, ...
    unitData.centroidY, unitData.nSpikes, unitData.firingRateHz, ...
    'VariableNames',{'unitId','memberElectrodes','centroidX','centroidY', ...
    'nSpikes','firingRateHz'});
exported.organoidLabel = unitData.organoidLabel;
exported.passesPairSpikeFilter = unitData.passesPairSpikeFilter;
end

function exported = roiExportTable(organoidRois)
% Row order carries no identity, so the type column is what says which
% organoid is which when this file is read back.
exported = table(string(organoidRois.type), organoidRois.label, ...
    organoidRois.position(:,1), organoidRois.position(:,2), ...
    organoidRois.position(:,3), organoidRois.position(:,4), ...
    'VariableNames', {'type','label','x','y','width','height'});
end
