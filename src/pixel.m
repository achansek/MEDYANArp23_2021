classdef pixel
%% This class containes information of pixels generated from MEDYAN frames.
    properties
        %Whether a pixel is empty or not
        occupiedstatus = false;
        %Has this pixel been accounted for in the breadth-first algorithm?
        accountedstatus = true;
        %number of monomers in the pixel
        count = 0;
        %mean concentration in muM of actin in the pixel
        meanconc = 0;
    end
    methods (Static)
        %Constructor
        function Npixels = setgetNpixels(data)
            persistent Var;
            if nargin
                Var = data;
            end
            Npixels = Var;
        end
        %Get the three dimensional coordinates of the pixel in the relative
        %axis based on 1 dimensional ID. For example, the first pixel with
        %center of mass [50, 50, 50] nm will have an idx1d of 1 and an
        %idx3d of [1 1 1].
        function idx3d = getidx3dfrom1d(idx1d)
            Npixelvec = pixel.setgetNpixels();
            crosssection = (Npixelvec(1)*Npixelvec(2));
            zaxispos = ceil(idx1d./crosssection);
            xyplanepos = idx1d - (zaxispos-1)*crosssection;
            yaxispos = ceil(xyplanepos/Npixelvec(1));
            xaxispos = xyplanepos - (yaxispos-1)*Npixelvec(1);
            idx3d = [xaxispos, yaxispos, zaxispos];
        end
        %Outputs the one dimensional IDs of the face touching neighbors of
        %a pixel.
        function neighborlist = getneighbors(idx1d)
            neighborlist=[];
            Npixels = pixel.setgetNpixels();
            crosssection = Npixels(1)*Npixels(2);
            %Get the three dimensional ID of pixel
            idx3d = pixel.getidx3dfrom1d(idx1d);
            %Move along X axis
            for i=-1:1
                %Reject if out of bounds
                if(idx3d(1)+i>Npixels(1) || idx3d(1)+i <= 0)
                    continue
                end
                %Move along Y axis
                for j =-1:1
                    %Reject if out of bounds
                    if(idx3d(2)+j>Npixels(2) || idx3d(2)+j <= 0)
                        continue
                    end
                    %Move along Z axis
                    for k=-1:1
                        %Reject if out of bounds
                        if(idx3d(3)+k>Npixels(3) || idx3d(3)+k <= 0)
                            continue
                        end
                        %Compute 3D ID of the pixel
                        neighbor3dindex = idx3d +[ i j k];
                        %Compute the 1D ID of the pixel
                        neighbor1dindex = (neighbor3dindex(3)-1)*crosssection + (neighbor3dindex(2)-1)*Npixels(1) + neighbor3dindex(1);
                        %Do not add self as a neighbor
                        if(neighbor1dindex~=idx1d)
                            neighborlist = [neighborlist, neighbor1dindex];
                        end
                    end
                end
            end
        end
    end
    methods
        %Alternate constructor
        function obj=pixel(N,varargin)
            if(nargin)
                for idx=1:N
                    obj(idx).occupiedstatus = false;
                end
            end
        end
        %Append pixel instance to an existing pixel array
        %obj - pixel array
        %obj2 - pixel instance to be inserted
        function obj=appendpixel(obj,obj2)
            if(nargin==2)
                obj(numel(obj)+1) = obj2;
                obj(numel(obj)).id = numel(obj);
            end
        end
    end
end