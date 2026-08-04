C65Target = 0.001; 
TimeVals = linspace(0, 48, 200);  % Test times between 0 and 48 hours
PlaceHolderTime = NaN;

for TimeInitGuess = TimeVals
    [C, ~] = NumericalFTCS_time(TimeInitGuess);
    % Check if concentration at C(65)/x=1000 reaches the target
    if C(65) >= C65Target
        PlaceHolderTime = TimeInitGuess;
        break;  % Exit the loop if the target is reached
    end
end

PlaceHolderTime

[C,y] = NumericalFTCS_time(PlaceHolderTime);
MinerPos750 = C(49);

[C,y] = NumericalFTCS_time(PlaceHolderTime+(250/360));
MinerPos1000 = C(65);

[C,y] = NumericalFTCS_time(PlaceHolderTime+(500/360));
MinerPos1250 = C(81);

[C,y] = NumericalFTCS_time(PlaceHolderTime+(750/360));
MinerPos1500 = C(97);

[C,y] = NumericalFTCS_time(PlaceHolderTime+(1000/360));
MinerPos1750 = C(113);

[C,y] = NumericalFTCS_time(PlaceHolderTime+(1250/360));
MinerPos2000 = C(129);