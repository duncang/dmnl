



Filename = 'data.txt';


filename_obsAMicro = fopen(Filename,'r');

TimeVec = [4, 5, 6,10,11,15, 22];
%for 15.00

k =  1;
j = 1;


Direction(1,1) = 0;
Speed(1,1) = 0;



Hgt = 5;  %5 ft



for i = 1:8000

    obslineMicro = fgetl(filename_obsAMicro);


    for m = 1:length(TimeVec)

        Time = TimeVec(m);

        if (str2num(obslineMicro(22:23)) == Time & str2num(obslineMicro(74:79)) == Hgt)

            if obslineMicro(53:55) == ' '

                Direction(m,k) = Direction(m,k-1);
            else

                Direction(m,k) = str2num(obslineMicro(53:55));
                
                if m ==length(TimeVec)
                k = k+1;
                end
            end

            if obslineMicro(46:49) == ' '

                Speed(m,j) = Speed(m,j-1);
            else

                Speed(m,j) = str2num(obslineMicro(46:49));
                                
                if m ==length(TimeVec)
                j = j+1;

                end
            end
        else    %if obslineMicro(
            %do nothing
            % Direction(i) =  Direction(i-1);
            %Speed(i) = Speed(i-1);
        end

        %end  %if (obslineMicro(22:23) ~= ' ' &

        %end %& obslineMicro(76:79) ~= ' ')
    end  % for m = 1:length(TimeVec)

end

fclose(filename_obsAMicro);