%% Function BAL_process.m
% Processes LTT balance and pressure data
% =========================================================================
% Tomas Sinnige - T.Sinnige@tudelft.nl 
% TU Delft - LR - AWEP - Flight Performance and Propulsion
%
% Version: 2.0
% Last updated:  17 Jan 2021
% First version: 17 October 2017
% =========================================================================
% | Version |    Date   |   Author  |              Changelog              |
% |---------|-----------|-----------|-------------------------------------|
% |   2.1   | 26/02/'21 |  P.Lopez  |  Removed PRS code                   |
% |---------|-----------|-----------|-------------------------------------|
% |   2.0   | 17/01/'21 |  P.Lopez  |  Change BAL and PRS to structure    |
% |         |           |    (TA)   |  type                               |
% |---------|-----------|-----------|-------------------------------------|
% |   1.0   | 28/05/'19 | T.Sinnige | Made code more generally applicable |
% |         |           |           | by removing part of previous code   |
% |         |           |           | related to isolated propeller setup.|
% |---------|-----------|-----------|-------------------------------------|
% |   0.0   | 17/10/'17 | T.Sinnige | First version                       |
% |---------|-----------|-----------|-------------------------------------|
% =========================================================================
% Inputs:  diskPath - root path on disk where data is stored
%          fnBAL    - filename(s) of the raw BAL data file(s)
%          fn0      - filename(s) of the zero measurement data file(s)
%          idxB     - structure containing indices balance data (W3D)
%          D        - propeller diameter [m]
%          S        - reference area [m^2]
%          b        - wing span [m]
%          c        - mean aerodynamic chord [m]
%          XmRefB   - moment reference points (x,y,z coordinates) in
%                     balance reference system [1/c] 
%          XmRefM   - moment reference points (x,y,z coordinates) in model 
%                     reference system [1/c] 
%          dAoA     - angle of attack offset (subtracted from measured 
%                     values)   
%          dAoS     - angle of sideslip offset (subtracted from measured 
%                     values)
%          modelType- type of wind-tunnel model
%          modelPos - orientation of wind-tunnel model
%          testSec  - LTT test-section number   
% -------------------------------------------------------------------------
% Outputs: BAL      - structure containing all the configurations.
%                     Each configuration contains an structure with the
%                     data of the different parameters.
% =========================================================================
function [BAL] = BAL_process(diskPath,fnBAL,fn0,idxB,D,S,b,c,XmRefB,XmRefM,dAoA,dAoS,modelType,modelPos,testSec,PRS)

%% Check inputs
% check whether diskPath ends with slash, and if not, append it
if ~strcmpi(diskPath(end),'/') || ~strcmpi(diskPath(end),'\')
    diskPath = [diskPath,'/'];    
end

% check whether inputs for raw data and zero data have the same length, and
% if not, throw an error
if length(fnBAL)~=length(fn0)
    error('Invalid inputs for fnBAL and fn0. Enter 1 zero offset file (fn0) for each raw data file (fnBAL).')
end

%% Process data
for i=1:length(fnBAL) % loop over the files that are to be loaded
    
    % extract configuration name
    BAL.config(i) = extractBetween(char(fnBAL{i}),'raw_','.txt');
    
    if isletter(BAL.config{i}(1))
        % first character is letter so can be used as field name -> continue
    else
        if strcmpi(BAL.config{i}(1),'+')
            BAL.config{i}(1) = 'p';
        elseif strcmpi(BAL.config{i}(1),'-')
            BAL.config{i}(1) = 'm';
        else
           error('Unexpected character used as first character of filename. Please add if statement here that catches this exception and changes the character into a letter.') 
        end
    end
    
    % remove plus/minus symbols from rest of config name (if applicable)  
    BAL.config{i} = strrep(BAL.config{i},'+','p');
    BAL.config{i} = strrep(BAL.config{i},'-','n');    

    % give status update
    display(['Processing balance data; configuration ''',BAL.config{i},'''; filename ''',fnBAL{i},'''.']);
    % load zero-measurement data (BAL0)
    BAL.windOff.(BAL.config{i}) = BAL_read0data(diskPath,fn0{i},idxB);
    % load measurement data (BAL)
    BAL.windOn.(BAL.config{i}) = BAL_readData([diskPath,fnBAL{i}],idxB);
    
    % check whether tunnel data were stored in the balance data file
    if all(BAL.windOn.(BAL.config{i}).pBar>0)
        % data included -> no problem
    else
        display('wind-tunnel data not stored in the balance data file -> get data from pressure data file.')
        
        tPRS = datetime(0,0,0,PRS.(PRS.config{i}).hr,PRS.(PRS.config{i}).min,PRS.(PRS.config{i}).sec);
        tBAL = datetime(0,0,0,BAL.windOn.(BAL.config{i}).hr,BAL.windOn.(BAL.config{i}).min,BAL.windOn.(BAL.config{i}).sec);

        % loop over each of the datapoints and get the correct data from
        % the PRS file
        for j=1:length(BAL.windOn.(BAL.config{i}).run)

            % check for corresponding measurement in PRS file
            idx = find(ismember(round(PRS.(PRS.config{i}).AoA*20)/20,round(BAL.windOn.(BAL.config{i}).AoA(j)*20)/20) & ... % angle of attack
                  ismember(round(PRS.(PRS.config{i}).AoS*20)/20,round(BAL.windOn.(BAL.config{i}).AoS(j)*20)/20) & ... % angle of sideslip
                  abs(BAL.windOn.(BAL.config{i}).rpsM1(j)-PRS.(PRS.config{i}).rpsM1) < 0.5 & ... % motor 1 rps 
                  abs(seconds(tBAL(j)- tPRS)) < 30); % datapoint time
              
            if length(idx)>1 && ismember(j,idx)
                [~,idxT]=min(abs(seconds(tBAL(j)- tPRS(idx)))) ;               
                idx = idx(idxT);
            end
            
            BAL.windOn.(BAL.config{i}).dPb(j) = PRS.(PRS.config{i}).dPb(idx);
 
        end

    end
        
    % compute calibrated balance forces and moments
    BAL.windOn.(BAL.config{i}) = BAL_forces(BAL.windOn.(BAL.config{i}),BAL.windOff.(BAL.config{i}),idxB,D,S,b,c,XmRefB,XmRefM,dAoA,dAoS,modelType,modelPos,testSec);    
    
end % end for loop over filenames

end % end of function BAL_process