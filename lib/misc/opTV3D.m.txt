classdef opTV3D < opSpot
    % opTV, anisotropic TV operator
    properties (Access = private)
        funHandle
        imSize
    end
    
    methods
        
        % constructor
        function op = opTV3D(imSize)
            
            if numel(imSize) == 1
                imSize = [imSize, imSize, imSize];
            elseif numel(imSize) ~= 3
                error('Invalid image dimensions!');
            end
            
            op = op@opSpot('opTV3D', prod(imSize-[1 0 0])+prod(imSize-[0 1 0])+prod(imSize-[0 0 1]), prod(imSize));
            
            op.funHandle = @opTV_internal;
            
            op.imSize = imSize;
            
            op.cflag = false;
            op.sweepflag = false;
        end
    end % Public methods
    
    methods(Access = protected)
        
        % Multiplication
        function y = multiply(op, x, mode)
            y = op.funHandle(op, x, mode);
        end
        
    end % Protected methods
    
    methods(Access = private)
        
        function y = opTV_internal(op, x, mode)
            
            % forward product
            if mode == 1
                
                % allocate output
                y = zeros(op.m,1);
                
                % reshape x to image dimensions
                x = reshape(x, op.imSize);
                
                % horizontal filtering
                temp = diff(x,1,1);
                y(1:prod(op.imSize-[1 0 0])) = temp;
                
                % vertical filtering
                temp = diff(x,1,2);
                y(prod(op.imSize-[1 0 0])+(1:prod(op.imSize-[0 1 0]))) = temp;
                
                % slice filtering
                temp = diff(x,1,3);
                y(prod(op.imSize-[1 0 0])+prod(op.imSize-[0 1 0])+1:end) = temp;
                
            % transposed product
            else
                x1 = reshape(x(1:prod(op.imSize-[1 0 0])),op.imSize-[1 0 0]);
                y1 = -diff(padarray(x1,[1 0 0],0,'both'),1,1);
                
                x2 = reshape(x(prod(op.imSize-[1 0 0])+(1:prod(op.imSize-[0 1 0]))),op.imSize-[0 1 0]);
                y2 = -diff(padarray(x2,[0 1 0],0,'both'),1,2);
                
                x3 = reshape(x(prod(op.imSize-[1 0 0])+prod(op.imSize-[0 1 0])+1:end),op.imSize-[0 0 1]);
                y3 = -diff(padarray(x3,[0 0 1],0,'both'),1,3);
                
                y = y1(:) + y2(:) + y3(:);
            end
            
        end
    end  
end
    