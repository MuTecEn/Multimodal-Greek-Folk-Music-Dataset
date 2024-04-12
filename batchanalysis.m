mirverbose(0)
close all

b=[];
filenames = {};

tab = readtable('keysignature.csv');

[filenames,patterns] = recurse(filenames,[],tab,'');
ons = {};

if 1
    %% Displaying the motifs in Command Window

    minlength = 3;

    for h = 1:length(filenames)
        filename = filenames{h};
        nmat = readmidi(filename);
        ons{h} = nmat(:,1);
    end

    f = 0;
    for i = 1:length(patterns)
        patt = patterns(i);

        found = 0;
        for j = 1:length(patt.occurrences)
            if isempty(patt.occurrences(j).extensions) && ...
                    patt.occurrences(j).pattern.length > minlength...
                    && isempty(patt.occurrences(j).cycle)
                found = 1;
                break
            end
        end
        if ~found
            continue
        end

        if patt.abstract
            continue
        end

        s.branches = patt;
        s.closedbranches = patt;
        s.general = patt;
        s.length = patt.length;
        if s.length > minlength
            coord = zeros(0,2);
            for h = 1:length(s.closedbranches)
                p = s.closedbranches(h);
                %                         if p.length < 4
                %                             continue
                %                         end
                occs = p.occurrences;
                firstfile = [];
                multi = 0;
                displayer = struct('fileID',{},'beat',{});
                for k = 1:length(occs)
                    fileID = occs(k).suffix.to.sequence.name;
                    if isa(occs(k).suffix,'pat.event')
                        N2 = occs(k).suffix.address;
                    else
                        N2 = occs(k).suffix.to.address;
                    end
                    nk = N2;

                    if isempty(N2)
                        continue
                    end

                    if isa(occs(k).suffix,'pat.event')
                        t2 = ons{fileID}(occs(k).suffix.address);
                    else
                        t2 = ons{fileID}(occs(k).suffix.to.address);
                    end
                    occk = occs(k);
                    patk = p;
                    while ~isempty(patk.parent) && ...
                            ~isempty(patk.parent.parent)
                        nk = [occk.suffix.to.address nk];
                        occk = occk.prefix;
                        patk = patk.parent;
                    end
                    N1 = occk.suffix;
                    if isa(N1,'pat.syntagm')
                        N1 = N1.to;
                    end
                    N1 = N1.address;
                    %                             if isempty(pool{end}.start)
                    %                                 pool{end}.start = N1;
                    %                             end

                    if isempty(N1)
                        continue
                    end

                    found = 0;
                    for l = 1:size(coord,1)
                        if N1 < coord(l,2) && N2 > coord(l,1)
                            found = 1;
                            N1 = min(N1,coord(l,1));
                            N2 = max(N2,coord(l,2));
                            break
                        end
                    end
                    if ~found
                        coord(end+1,:) = [N1,N2];
                    end

                    if isa(occk.suffix,'pat.event')
                        t1 = ons{fileID}(occk.suffix.address);
                    else
                        t1 = ons{fileID}(occk.suffix.to.address);
                    end

                    if t1 > t2
                        warning('Occurrence error');
                        break
                    end

                    if isempty(firstfile)
                        firstfile = fileID;
                    elseif fileID ~= firstfile
                        multi = 1;
                    end

                    displayer(end+1).fileID = fileID;
                    displayer(end).beat = t1;

                end

                if multi
                    f = f + 1;
                    fprintf(['== ',num2str(f),'\n']);
                    desc = p.display; %(0);
                    for k = 1:length(displayer)
                        disp(['Tune',filenames{displayer(k).fileID},', Beat ',num2str(displayer(k).beat)])
                    end
                end
            end
        end
    end
else
    %% Displaying the motifs

    minlength = 2;

    for h = 1:length(filenames)
        figure
        hold on

        filename = filenames{h};
        nmat = readmidi(filename);
        ons{h} = nmat(:,1);
        chro = nmat(:,4);

        for i = 1:length(ons)
            plot([ons{h}(i),ons{h}(i)],[chro(i),chro(i)],'Color','k','LineWidth',1.5);
            plot(ons{h}(i),chro(i),'dk',...
                'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7);
        end

        xl = xlim;

        y = min(chro);
        stafflines = [43 47 50 53 57 64 67 71 74 77];
        stafflines(stafflines < y) = [];

        for i = stafflines
            plot(xl,[i i],'Color',[.7 .7 .7],'LineWidth',2)
        end

        for i = 0:4:xl(end)
            plot([i i],ylim,'Color',[.7 .7 .7])
        end

        title(filename)
    end

    f = 0;
    y = y - 1;
    for i = 1:length(patterns)
        patt = patterns(i);

        found = 0;
        for j = 1:length(patt.occurrences)
            if isempty(patt.occurrences(j).extensions) && ...
                    patt.occurrences(j).pattern.length > minlength...
                    && isempty(patt.occurrences(j).cycle)
                found = 1;
                break
            end
        end
        if ~found
            continue
        end

        if patt.abstract
            continue
        end

        f = f + 1;

        s.branches = patt;
        s.closedbranches = patt;
        s.general = patt;
        s.length = patt.length;
        if s.length > minlength
            coord = zeros(0,2);
            fprintf(['== ',num2str(f),'\n']);
            if f > 1
                y = y - .3;
                %                         line([0 10],[y y],'Color','k','LineWidth',2);
                %                         y = y - .5;
            end
            text(0,y+.3,num2str(f),'Color','k','VerticalAlignment','Top');

            for h = 1:length(s.closedbranches)
                p = s.closedbranches(h);
                %                         if p.length < 4
                %                             continue
                %                         end
                if h > 1
                    y = y-.3;
                end
                desc = p.display; %(0);
                %fprintf([desc,'\n']);

                occs = p.occurrences;
                for k = 1:length(p.specific)
                    for k2 = 1:length(p.specific(k).occurrences)
                        found = 0;
                        for k3 = 1:length(occs)
                            if isequal(occs(k3).suffix,...
                                    p.specific(k).occurrences(k2))
                                found = 1;
                            end
                        end
                        if ~found
                            if isempty(occs)
                                occs = p.specific(k)...
                                    .occurrences(k2);
                            else
                                occs(end+1) = p.specific(k)...
                                    .occurrences(k2);
                            end
                        end
                    end
                end

                col = num2col(f);

                for k = 1:length(occs)
                    fileID = occs(k).suffix.to.sequence.name;

                    figure(fileID)

                    if isa(occs(k).suffix,'pat.event')
                        N2 = occs(k).suffix.address;
                    else
                        N2 = occs(k).suffix.to.address;
                    end
                    nk = N2;

                    if isempty(N2)
                        continue
                    end

                    if isa(occs(k).suffix,'pat.event')
                        t2 = ons{fileID}(occs(k).suffix.address);
                    else
                        t2 = ons{fileID}(occs(k).suffix.to.address);
                    end
                    occk = occs(k);
                    patk = p;
                    while ~isempty(patk.parent) && ...
                            ~isempty(patk.parent.parent)
                        nk = [occk.suffix.to.address nk];
                        occk = occk.prefix;
                        patk = patk.parent;
                    end
                    N1 = occk.suffix;
                    if isa(N1,'pat.syntagm')
                        N1 = N1.to;
                    end
                    N1 = N1.address;
                    %                             if isempty(pool{end}.start)
                    %                                 pool{end}.start = N1;
                    %                             end

                    if isempty(N1)
                        continue
                    end

                    found = 0;
                    for l = 1:size(coord,1)
                        if N1 < coord(l,2) && N2 > coord(l,1)
                            found = 1;
                            N1 = min(N1,coord(l,1));
                            N2 = max(N2,coord(l,2));
                            break
                        end
                    end
                    if ~found
                        coord(end+1,:) = [N1,N2];
                    end

                    if isa(occk.suffix,'pat.event')
                        t1 = ons{fileID}(occk.suffix.address);
                    else
                        t1 = ons{fileID}(occk.suffix.to.address);
                    end

                    if t1 > t2
                        warning('Occurrence error');
                        break
                    end
                    line([t1 t2],[y y],'Marker','d','Color',col,...
                        'LineWidth',2,'MarkerSize',6);

                end
            end
        end
    end

end








function [filenames,patterns] = recurse(filenames,patterns,tab,path)
    d = dir;
    for i = 3:5 %length(d)
        if d(i).isdir
            disp('/////////////////')
            d(i).name
            cd(d(i).name);
            if isempty(path)
                path2 = d(i).name;
            else
                path2 = [path,'/',d(i).name];
            end
            [filenames,patterns] = recurse(filenames,patterns,tab,path2);
            cd ..
        else
            d(i).name
            signature = [];
            for j = 1:size(tab,1)
                if strcmp(tab{j,1}{1},d(i).name)
                    signature = tab{j,2:6};
                    beatsinbar = tab{j,7};
                    break
                end
            end
            if isempty(signature)
                warning('No signature found!!')
                signature = [NaN NaN NaN NaN NaN];
                beatsinbar = 4;
            end
            if ~isempty(signature)
                [patterns,ok] = analyzeMIDI(d(i).name,signature,beatsinbar,patterns,length(filenames)+1);
                if ok
                    if isempty(path)
                        filenames{end+1} = [d(i).name];
                    else
                        filenames{end+1} = [path,'/',d(i).name];
                    end
                end
            end
        end
    end
end





