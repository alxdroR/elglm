function options = parseNameValueoptions(options,varargin)
%parseNameValueoption - return structure of user-controlled options
% options = parseNameValueoptions(options,varargin)
%
%   Use this file to parse Name/Value pairs when a function has
% optional input.  If directly sending varargin from funciton remember
% to call parseNameValueoption as  parseNameValueoptions(options,varargin{:})
%
% Ex: you have a function called myfun(x,varargin) that allows for the optional inputs
%     called op1 and op2. someone makes the call
%
%     myfun(x,'op1',3)
%
%    Inside myfun you must have this
%    options = {'op1',[],'op1',[]}
%    or default options instead of empty arrays
%    Inside my fun add
%    options = parseNameValueoptions(options,varargin{:});
%   Now options is a structure with field op1 set to 3 and
%   and op2 set to the default [];
%
% read the acceptable names
optionNames = fieldnames(options);

%# count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
    error('optional inputs need  propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
    inpName = pair{1};
    %inpName = lower(pair{1}); %# make case insensitive
    
    if any(strcmp(inpName,optionNames))
        %# test for the right class here
        options.(inpName) = pair{2};
    else
        %   error('%s is not a recognized parameter name',inpName)
    end
end
end

