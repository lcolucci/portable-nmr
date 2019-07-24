%% TextLocation
% To automatically place a text object at the 'Best' location similar to a 'legend' object, use the attached MATLAB function 'TextLocation'. To use this, execute the following steps:
% 
% 1. Download the attached function and store it in a location that is in MATLAB path.
% 
% 2. Call this function with required text after creating the figure. For example:
% 
% a. Insert the text 'Hello' in the 'best' location in the current figure:
% 
% plot(1:10);
% TextLocation('Hello','Location','best');
% b. Insert the text 'World' in the 'southwest' location in the current figure:
% 
% TextLocation('World','Location','southwest');
% Please note that the use of this function will delete any legend currently in the figure window. If you wish to use a legend in your figure, please specify it after calling the 'TextLocation' function.
%
% Taken from: https://www.mathworks.com/matlabcentral/answers/98082-how-can-i-automatically-specify-a-best-location-property-for-the-textbox-annotation-function

function hOut = TextLocation(textString,varargin)

l = legend(textString,varargin{:});
t = annotation('textbox');
t.String = textString;
t.Position = l.Position;
t.FontSize = l.FontSize; 
delete(l);
t.LineStyle = 'None';

if nargout
    hOut = t;
end
end