function tests = testConnectionRuleSettings
% The launcher's settings block is the connection rule. Every setting in it
% has to reach the analysis, and the values have to be the calibrated ones:
% a setting that is never passed is decoration, and a value that drifts makes
% one recording's connection counts incomparable with the next.
tests = functiontests(localfunctions);
end


% Every setting that decides which pairs become connections.
function names = ruleSettings()
names = ["MinElectrodeFiringRateHz", "MaxElectrodeFiringRateHz", ...
         "DbscanRadiusUm", "DbscanMinElectrodes", "MergeWindowMs", ...
         "MinUnitSpikesForPair", "CorrelationBinMs", "LatencyRangeMs", ...
         "PeakLatencyRangeMs", "MinPeakCount", "PeakToBaselineMultiplier", ...
         "MaxPeakEfficacy", "DuplicateRadiusUm", "MaxNearZeroLagRatio"];
end


function setupOnce(testCase)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
f = fullfile(root, "connectoid_launcher.m");
assert(isfile(f), 'No launcher at %s', f);
testCase.TestData.launcher = string(fileread(f));
end


function testEverySettingReachesTheAnalysis(testCase)
text = testCase.TestData.launcher;
for s = ruleSettings()
    testCase.verifyNotEmpty(settingValue(text, s), ...
        sprintf('connectoid_launcher does not set config.%s', s));
    % The call has to hand it over; setting it is not applying it.
    passed = ~isempty(regexp(text, "'" + s + "'\s*,\s*config\." + s, 'once'));
    testCase.verifyTrue(passed, sprintf( ...
        'config.%s is set but never passed to the analysis, so editing it does nothing', s));
end
end


function testTheValuesAreTheOnesInTheMethods(testCase)
% Changing one of these means re-running every recording in the dataset.
settled = struct('MinElectrodeFiringRateHz', "0.01", 'MinUnitSpikesForPair', "50", ...
                 'MinPeakCount', "10", 'PeakToBaselineMultiplier', "5", ...
                 'MaxPeakEfficacy', "0.5", 'DuplicateRadiusUm', "60", ...
                 'MaxNearZeroLagRatio', "2", 'DbscanRadiusUm', "60", ...
                 'MergeWindowMs', "1.5", 'CorrelationBinMs', "0.5");
for s = string(fieldnames(settled))'
    testCase.verifyEqual(settingValue(testCase.TestData.launcher, s), ...
        settled.(char(s)), sprintf('config.%s is not the settled value', s));
end
end


function v = settingValue(text, name)
% The right-hand side of "config.<name> = ...;", trailing comment removed.
m = regexp(text, "(?m)^config\." + name + "\s*=\s*([^;]*);", 'tokens', 'once');
if isempty(m)
    v = "";
else
    v = strtrim(string(m{1}));
end
end
