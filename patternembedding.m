minlength = 3;
s = [];

nfiles = length(filenames);
simatrix = zeros(nfiles);
lengths = zeros(1,nfiles);

for h = 1:nfiles
    nmat = readmidi(filenames{h});
    lengths(h) = size(nmat,1);
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
            counts = zeros(1,length(filenames));
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

                counts(fileID) = counts(fileID) + 1;
            end

            plength = patt.occurrences(1).length;
            active = find(counts);
            for j = 1:length(active)
                for k = j+1:length(active)
                    file1 = active(j);
                    file2 = active(k);
                    weight1 = counts(file1) * plength / lengths(file1);
                    weight2 = counts(file2) * plength / lengths(file2);
                    simatrix(file1,file2) = simatrix(file1,file2) ...
                        + weight1 * weight2;
                    simatrix(file2,file1) = simatrix(file1,file2);
                end
            end
        end
    end
end

dismatrix = 1./simatrix;
for i = 1:size(dismatrix,1)
    dismatrix(i,i) = 0;
end
dismatrix(isinf(dismatrix)) = 10000;

figure
if 0
    y = mytsne(dismatrix,'Algorithm','Exact');
    s = scatter(y(:,1),y(:,2),'o','filled');
    s.SizeData = 100;
    for i = 1:size(y,1)
        n = filenames{i};
        [filepath,name,ext] = fileparts(n);
        n = name;

        text(y(i,1),y(i,2),n,'FontSize',15,'Color',[.6 .6 .6])
    end
else
    y = mytsne(dismatrix,'Algorithm','Exact','NumDimensions', 3);
    s = scatter3(y(:,1),y(:,2),y(:,3),'o','filled');
    s.SizeData = 100;
    for i = 1:size(y,1)
        n = filenames{i};
        [filepath,name,ext] = fileparts(n);
        n = name;

        text(y(i,1),y(i,2),y(i,3),n,'FontSize',15,'Color',[.6 .6 .6])
    end
end