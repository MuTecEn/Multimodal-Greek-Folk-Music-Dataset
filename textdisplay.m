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