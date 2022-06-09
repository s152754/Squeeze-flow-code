clear all; close all; clc;
% Load data
iframes = importdata('InitialFrames.txt');
icase           = 1; %F = 1, F025 = 2, F05 = 3, V02 = 4
SecFrame        = 4;
Frames_Between  = SecFrame;
Frame_init      = iframes.data(:,icase);
Frame_middle    = Frame_init+Frames_Between;

crop1           = 1; 
crop2           = 1080; 
crop3           = 250; 
crop4           = 1400;

sens = 0.9;
ID = 1;
for i = 1:3 %1:length(Frame_init)
    filename        = ["videos/GlycerolT_V01_" + num2str(i) + ".mp4"] %#ok<NOPTS> 
    v               = VideoReader(filename); %#ok<TNMLP> 

    readFramei      = read(v,Frame_init(i)); 
    Frames_Between = SecFrame;
    ID = 1;
    while Frames_Between < v.NumFrames-Frame_init(i)
        Frames_Between
        Frame_middle    = Frame_init+Frames_Between;
        readFramem      = read(v,Frame_middle(i));
        Frame_i         = readFramei(crop1:crop2,crop3:crop4);
        Frame_m         = readFramem(crop1:crop2,crop3:crop4);    
        Frame_sub       = imsubtract(Frame_i,Frame_m);
        [BW,intensity]  = createBWimage(Frame_sub);

        if Frames_Between == SecFrame
            [CentIn, Rin] = imfindcircles(BW,[100 150],'ObjectPolarity','dark','Sensitivity',0.92);
        end

        if Frames_Between < 20
            Rmin = 100;  Rmax = 150;
        elseif Frames_Between > 20 && Frames_Between< 500
            Rmin = 150; Rmax = 300;
        elseif Frames_Between < 5000
            Rmin = 200; Rmax = 500;
        elseif Frames_Between < 20000
            Rmin = 300; Rmax = 600;
        elseif Frames_Between < v.NumFrames
            Rmin = 300; Rmax = 600;
        end

        [centersBright, radiiBright] = imfindcircles(BW,[Rmin Rmax],'ObjectPolarity','dark','Sensitivity',sens);

        while isempty(centersBright)
            sens  = sens + 0.01;
            [centersBright, radiiBright] = imfindcircles(BW,[Rmin Rmax],'ObjectPolarity','dark','Sensitivity',sens);
        end  
        
%         figure()
%         imshow(BW)
%         hold on
%         viscircles(centersBright, radiiBright,'Color','b');

        [muptl,sigptl] = pixeltolength(); 

        Rout(i,ID)    = radiiBright*muptl;
        CentOut(i,ID,1:2) = centersBright*muptl;
        time(i,ID) = Frames_Between/v.FrameRate;
        
        Frames_Between = Frames_Between*2;
        sens = 0.9;
        ID = ID + 1;        
    end 
    Rin         = Rin*muptl;
    CentIn      = CentIn*muptl;

    
    Radius(i,:)  = [Rin Rout(i,:)];
    Time(i,:)    = [0 time(i,:)];

end




%% check with old measurements
Rold = importdata("Gly_V01_F.txt");
muT = Rold.data(:,1);
muR = Rold.data(:,2);
stdR = Rold.data(:,3);

figure()
plot(Time(1,:),Radius(1,:),'Color','-ob')
hold on
plot(Time(2,:),Radius(2,:),'Color','-og')
plot(Time(3,:),Radius(3,:),'Color','-or')
errorbar(muT,muR,2*stdR,'k')



