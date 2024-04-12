function [patterns,ok] = analyzeMIDI(filename,spell,beatsinbar,patterns, fileID)

if nargin < 3
    beatsinbar = 4;
end
if nargin < 4
    patterns = [];
    fileID = 0;
end

pat.verbose(0);

%% Loading the MIDI file
try
    nmat = readmidi(filename);
    ok = 1;
catch
    ok = 0;
    return
end

chro = nmat(:,4);
chro7 = mod(chro,12);
oct = floor(chro/12);
dia = chro7;
if strcmp(filename,'Syrtos_Nysyros.mid')
    for i = 1:length(chro7)
        switch chro7(i)
            case 0
                dia(i) = 6;
            case 1
                dia(i) = 0;
            case 2
                dia(i) = 1;
            case 3
                dia(i) = 1;
            case 4
                dia(i) = 1;
            case 5
                dia(i) = 2;
            case 6
                dia(i) = 3;
            case 7
                dia(i) = 3;
            case 8
                dia(i) = 4;
            case 9
                dia(i) = 4;
            case 10
                dia(i) = 5;
            case 11
                dia(i) = 6;
        end
        dia(i) = dia(i) + 7*oct(i);
    end
else
    for i = 1:length(chro7)
        switch chro7(i)
            case 0
                dia(i) = 0;
            case 1
                dia(i) = spell(1);
            case 2
                dia(i) = 1;
            case 3
                dia(i) = spell(2);
            case 4
                dia(i) = 2;
            case 5
                dia(i) = 3;
            case 6
                dia(i) = spell(3);
            case 7
                dia(i) = 4;
            case 8
                dia(i) = spell(4);
            case 9
                dia(i) = 5;
            case 10
                dia(i) = spell(5);
            case 11
                dia(i) = 6;
        end
        dia(i) = dia(i) + 7*oct(i);
    end
end

ons = nmat(:,1);
beat = mod(ons,beatsinbar);
bar = floor(ons/beatsinbar) + 1;
beat(beat > 3.8) = 0; % due to rounding error issues
beat = round(beat * beatsinbar) / beatsinbar;

notes = struct('chro',num2cell(chro),'dia',num2cell(dia),...
    'ons',num2cell(ons),'beat',num2cell(beat),'bar',num2cell(bar),...
    'dur',0,'chan',num2cell(nmat(:,3)),'filename',filename,'fileID',fileID);

chro = [notes.chro];
ons = [notes.ons];
dur = [notes.dur];
off = ons + dur;

ps = seq.paramstruct('MIDI',{'dia','beat'},1,@validfield);

paramdia = seq.paramtype('dia');
paramdia.inter = seq.paraminter(@(x,y) x-y);
ps = ps.setfield('dia',paramdia);

parambeat = seq.paramtype('beat');
ps = ps.setfield('beat',parambeat);

%% Actual pattern analysis
[sequence, patterns] = analyze(notes,ps,@parseSymbol,patterns);

if nargin < 4
    %% Displaying the motifs

    minlength = 2;

    figure
    hold on

    for i = 1:length(notes)
        plot([ons(i),off(i)],[chro(i),chro(i)],'Color','k','LineWidth',1.5);
        plot(ons(i),chro(i),'dk',...
            'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',7);
    end

    xl = xlim;

    y = min(chro);
    stafflines = [43 47 50 53 57 64 67 71 74 77];
    stafflines(stafflines < y) = [];

    for i = stafflines
        plot(xl,[i i],'Color',[.7 .7 .7],'LineWidth',2)
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

%         if inter
%             thissong = 0;
%             othersongs = 0;
%             for j = 1:length(patt.occurrences)
%                 namej = patt.occurrences(j).suffix.to.sequence.name;
%                 if strcmp(namej,filename)
%                     thissong = 1;
%                 else
%                     othersongs = 1;
%                 end
%                 if thissong && othersongs
%                     break
%                 end
%             end
%             if ~thissong || ~othersongs
%                 continue
%             end
%         end
% 
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
%                 for k = 1:length(p.specific)
%                     for k2 = 1:length(p.specific(k).occurrences)
%                         found = 0;
%                         for k3 = 1:length(occs)
%                             if isequal(occs(k3).suffix,...
%                                     p.specific(k).occurrences(k2))
%                                 found = 1;
%                             end
%                         end
%                         if ~found
%                             if isempty(occs)
%                                 occs = p.specific(k)...
%                                     .occurrences(k2);
%                             else
%                                 occs(end+1) = p.specific(k)...
%                                     .occurrences(k2);
%                             end
%                         end
%                     end
%                 end

                col = num2col(f);

                for k = 1:length(occs)
%                     if ~strcmp(occs(k).suffix.to.sequence.name,filename)
%                         continue
%                     end

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
                        t2 = ons(occs(k).suffix.address);
                    else
                        t2 = ons(occs(k).suffix.to.address);
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
                        t1 = ons(occk.suffix.address);
                    else
                        t1 = ons(occk.suffix.to.address);
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

    for i = 0:4:xl(end)
        plot([i i],ylim,'Color',[.7 .7 .7])
    end

    drawnow

end


function [p,txt] = parseSymbol(note,ps)
% The parametrical description of the new event
p = ps.type2val; % Defined from the parametrical space ps
% ... by associating a value to each parameter (here just one, the
% "dimension")
p = p.setfield('dia',seq.paramval(ps.getfield('dia'),note.dia));
p = p.setfield('beat',seq.paramval(ps.getfield('beat'),note.beat))
txt = note.chro;

function y = validfield(i,options)
y = 1;