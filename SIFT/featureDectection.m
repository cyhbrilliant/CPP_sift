function varargout = featureDectection(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @featureDectection_OpeningFcn, ...
                   'gui_OutputFcn',  @featureDectection_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end


function featureDectection_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;


guidata(hObject, handles);




function varargout = featureDectection_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


function pushbutton_Load_Callback(hObject, eventdata, handles)

filenameA = 'test10.jpg';
filenameB = 'test11.jpg';
imgA = imread(filenameA);
imgB = imread(filenameB);
ImageInformationA = imfinfo(filenameA);
if strcmp(ImageInformationA.ColorType,'truecolor')
    hcsc = vision.ColorSpaceConverter('Conversion','RGB to intensity');
    imgA = step(hcsc, imgA);
end
ImageInformationB = imfinfo(filenameB);
if strcmp(ImageInformationB.ColorType,'truecolor')
    hcsc = vision.ColorSpaceConverter('Conversion','RGB to intensity');
    imgB = step(hcsc, imgB);
end
imshowpair(imgA,imgB,'montage');
handles.imgA = imgA;
handles.imgB = imgB;
guidata(hObject,handles);


function pushbutton_detect_Callback(hObject, eventdata, handles)

feature = get(handles.popupmenu_feature,'Value');
imgA = handles.imgA;
imgB = handles.imgB;
switch feature
    case 2           %'FAST'
        pointsA = detectFASTFeatures(imgA);
        pointsB = detectFASTFeatures(imgB);
    case 3           %'Harris'
        pointsA = detectHarrisFeatures(imgA);
        pointsB = detectHarrisFeatures(imgB);
    case 4          %'MinEigen'
        pointsA = detectMinEigenFeatures(imgA);
        pointsB = detectMinEigenFeatures(imgB);
    case 5           %'MSER'
        pointsA = detectMSERFeatures(imgA);
        pointsB = detectMSERFeatures(imgB);
    case 6           %'SURF'
        pointsA = detectSURFFeatures(imgA);
        pointsB = detectSURFFeatures(imgB);
    otherwise
        pointsA = detectSURFFeatures(imgA);
        pointsB = detectSURFFeatures(imgB);
end
locA = pointsA.Location;
locB = pointsB.Location;
imshowpair(imgA,imgB,'montage');
hold on
for i = 1: length(locA)
    plot(locA(i,1),locA(i,2),'*')
end
for i = 1: length(locB)
    plot(locB(i,1)+size(imgB,2),locB(i,2),'ro')
end
hold off
handles.pointsA = pointsA;
handles.pointsB = pointsB;
handles.feature = feature;
guidata(hObject, handles);


function pushbutton_extract_Callback(hObject, eventdata, handles)

pointsA = handles.pointsA;
pointsB = handles.pointsB;
imgA = handles.imgA;
imgB = handles.imgB;
switch handles.feature
    case 2           %'FAST'
        [fA, vptsA] = extractFeatures(imgA, pointsA);
        [fB, vptsB] = extractFeatures(imgB, pointsB);
    case 3           %'Harris'
        [fA, vptsA] = extractFeatures(imgA, pointsA);
        [fB, vptsB] = extractFeatures(imgB, pointsB);
    case 4          %'MinEigen'
        [fA, vptsA] = extractFeatures(imgA, pointsA);
        [fB, vptsB] = extractFeatures(imgB, pointsB);
    case 5           %'MSER'
        [fA, vptsA] = extractFeatures(imgA, pointsA);
        [fB, vptsB] = extractFeatures(imgB, pointsB);
    case 6           %'SURF'
        [fA, vptsA] = extractFeatures(imgA, pointsA);
        [fB, vptsB] = extractFeatures(imgB, pointsB);
    otherwise
        [fA, vptsA] = extractFeatures(imgA, pointsA);
        [fB, vptsB] = extractFeatures(imgB, pointsB);
end
handles.fA = fA;
handles.fB = fB;
handles.vptsA = vptsA;
handles.vptsB = vptsB;
guidata(hObject, handles);


function pushbutton_save_Callback(hObject, eventdata, handles)




function pushbutton_match_Callback(hObject, eventdata, handles)

if isempty(handles.fA)
    return;
end
fA = handles.fA;
fB = handles.fB;
vptsA = handles.vptsA;
vptsB = handles.vptsB;
imgA = handles.imgA;
imgB = handles.imgB;

method = get(handles.popupmenu_matchMethod,'Value');
switch method
    case 2
        indexPairs = matchFeatures(fA, fB,'Method','Threshold', 'Prenormalized', false);
    case 3
        indexPairs = matchFeatures(fA, fB,'Method','NearestNeighborSymmetric', 'Prenormalized', false);
    case 4
        indexPairs = matchFeatures(fA, fB,'Method','NearestNeighborRatio', 'Prenormalized', false);
    otherwise
        indexPairs = matchFeatures(fA, fB,'Method','NearestNeighborRatio', 'Prenormalized', false);
end
        

matched_ptsA = vptsA(indexPairs(:,1));
matched_ptsB = vptsB(indexPairs(:,2));
matched_locA = matched_ptsA.Location;
matched_locB = matched_ptsB.Location;
matched_locB(:,1) = matched_locB(:,1) + size(imgB,2);
% Show matched points
imshowpair(imgA,imgB,'montage');
hold on
plot(matched_locA(:,1),matched_locA(:,2),'*',matched_locB(:,1),matched_locB(:,2),'ro')
for i = 1: size(matched_locA,1)
    line([matched_locA(i,1),matched_locB(i,1)],[matched_locA(i,2),matched_locB(i,2)],...
        'Color','g');
end
hold off

function pushbutton_close_Callback(hObject, eventdata, handles)

close

% --- Executes on selection change in popupmenu_feature.
function popupmenu_feature_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function popupmenu_feature_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_matchMethod.
function popupmenu_matchMethod_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_matchMethod_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
