function varargout = NormaliztionExt(varargin)
nbIN     = nargin;
stdINPUT = true;
if nbIN>0 , stdINPUT = ~ischar(varargin{1}); end
if stdINPUT
    narginchk(6,7);
    X1 = varargin{1}; 
    X2 = varargin{2};
    wname = varargin{3};
    level = varargin{4};
    AfusMeth = varargin{5};
    DfusMeth = varargin{6};
    if nargin>6
        if strcmpi(varargin{7},'plot')
            flagPlot = true; 
        else
            flagPlot = false;
        end
    else
            flagPlot = false;
    end
else
    wname_DEF    = 'db1';
    level_DEF    = 2;
    fusMeth_DEF  = struct('name','linear','param',0.5);
    flagPlot_DEF = true;
    wname = wname_DEF;
    level = level_DEF;
    AfusMeth = fusMeth_DEF;
    DfusMeth = fusMeth_DEF;
    flagPlot = flagPlot_DEF;
    k    = 1;
    while k<=nbIN
        switch varargin{k}
            case 'X1'       , X1       = varargin{k+1};
            case 'X2'       , X2       = varargin{k+1};
            case 'wname'    , wname    = varargin{k+1};
            case 'level'    , level    = varargin{k+1};
            case 'AfusMeth' , AfusMeth = varargin{k+1};
            case 'DfusMeth' , DfusMeth = varargin{k+1};    
            case 'flagPlot' , flagPlot = varargin{k+1};
            otherwise
                error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'));
        end
        k = k+2;
    end
    if isempty(X1) || isempty(X2)
        error(message('Wavelet:FunctionArgVal:Invalid_X1X2Val'));
    end
    if isempty(wname)    , wname = wname_DEF;  end
    if isempty(level)    , level = level_DEF;  end
    if isempty(AfusMeth)  , AfusMeth = fusMeth_DEF;  end
    if isempty(DfusMeth)  , DfusMeth = fusMeth_DEF;  end
    if strcmp(flagPlot,'noplot') || isempty(flagPlot)
        flagPlot = false;
    end
    if ischar(X1) , dummy = load(X1); X1 = dummy.X; end
    if ischar(X2) , dummy = load(X2); X2 = dummy.X; end
end
if sum(size(X1)==size(X2))<2
    error(message('Wavelet:FunctionArgVal:Invalid_ImgSiz'));
end
tIMG1 = wfustree(X1,level,wname);
clear X1
tIMG2 = wfustree(X2,level,wname);
clear X2
[XFus,tFus] = wfusdec(tIMG1,tIMG2,AfusMeth,DfusMeth);
if flagPlot
    plot(tIMG1); plot(tIMG2); plot(tFus);
end
switch nargout
    case 0 ,
    case {1,2} , varargout = {XFus , tFus};
    otherwise  , varargout = {XFus , tFus, tIMG1,tIMG2};
end
