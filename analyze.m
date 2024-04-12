% Compact version of pat.analyze

function [outSequence,patterns] = analyze(inSequence,ps,parseSymbol,patterns)
if nargin < 3
    parseSymbol = @defaultParse;
end

options.explicit = true;
pat.cyclic(0);

if nargin < 2
    % Definition of the parametrical space. Here just one dimension.
    ps = seq.paramstruct('mystruct',{'dimension'},1);
    dimension = seq.paramtype('dimension');
    ps = ps.setfield('dimension',dimension);
end

if nargin < 4 || isempty(patterns)
    % Creating the root of the pattern tree
    root = pat.pattern([],[],[],ps);

    % Creating the "occurrence" of the root. All pattern occurrences will
    % extend from that "root occurrence".
    occ0 = root.occurrence([],[]);

    patterns = root;
else
    root = patterns(1);
    occ0 = root.occurrences(1);
end

% % Creating the note pattern: each note is an occurrence of the note
% % pattern. It is used to find interval patterns, since the first note has
% % no interval description.
% notepattern = pat.pattern([],root,ps);
% notepattern.general = root;
% root.specific = notepattern;

% The sequence progressively created, starts empty
outSequence = seq.Sequence([],inSequence(1).fileID);

previous = [];
currentbar = 0;
for i = 1:length(inSequence)
    % Each event successively considered in the sequence
    if inSequence(i).bar > currentbar
        disp('-------------')
        currentbar = inSequence(i).bar
    end

    [p,txt] = parseSymbol(inSequence(i),ps);
%     [i,txt]

    % The event is then instantiated as a pat.event object
    event = pat.event(outSequence,p);
    event.address = i;

    % ... integrated into the sequence under creation.
    outSequence = outSequence.integrate(event);

    % Now the two core mechanisms for pattern analysis:
    if ~isempty(previous)
        % 1. Trying to extent any pattern occurrence ending at the
        % previous event

        succ = pat.syntagm(previous,event,root,0);

        occ = previous.occurrences;
        for j = 1:length(occ)
            % Memorise the extension of the pattern with this new event
            patterns = memorize(occ(j),succ,root,options,patterns);
        end
    end

    % 2. Checking if the new event can be the start of a new pattern
    % occurrence (i.e., by extending the "root occurrence" occ0)
    patterns = memorize(occ0,event,root,options,patterns);
    % And here also, it consists in calling pat.occurrence.memorize, where
    % the occurrence to extend is the "root occurrence".

    previous = event;
end

%%
function patterns = memorize(occ,succ,root,options,patterns)
objpat = occ.pattern;

% Pattern recognition
% Can the pattern occurrence occ(j) be extended into an
% occurrence of a known child of the pattern?
recognize(objpat,occ,succ,root,options);

% Pattern discovery
% Can the pattern occurrence occ(j) be extended into an
% occurrence of a new child of the pattern?

% Corresponding to objpat.memory.fields{1}.learn(succ.parameter.fields{1},occ,succ,objpat,[],0,root);
param = succ.parameter;
% if length(param.fields) > 1
%     param.fields{2}.value = [];
% end

if isa(succ,'pat.event')
    occurrences = succ.occurrences;
else
    occurrences = succ.to.occurrences;
end

    patterns = learn(objpat.memory.fields{1}.inter,param.fields{1}.inter,param,occ,succ,objpat,root,true,patterns);
    patterns = learn(objpat.memory.fields{2},param.fields{2},param,occ,succ,objpat,root,true,patterns);


%%
function recognize(obj,occ,succ,root,options)
if isempty(obj.parent)
    addr = 0;
else
    addr = obj.address;
end
if ismember(addr,succ.processed)
    return
end
if isempty(succ.processed)
    succ.processed = addr;
else
    succ.processed(end+1) = addr;
end

% We look at each child
for i = 1:length(obj.children)
    child = obj.children{i};
    param = child.parameter;

    stop = false;
    if isa(succ,'pat.syntagm')
        to = succ.to;
    else
        to = succ;
    end

    % Generalizing if necessary.
    if isempty(occ.suffix)
        param0 = [];
    else
        param0 = occ.suffix.parameter;
    end
    [test, common] = succ.parameter.implies(param,param0);

    if ~test
        if isempty(common.fields{2}.value) || ...
                (~isequal(obj,root) && ...
                isempty(common.fields{1}.inter.value))
            continue
        end
        param = common;
    end

    if ~isempty(occ)
        % Don't go further if extension already exists.
        found = 0;
        for j = 1:length(occ.extensions)
            if occ.extensions(j).parameter.implies(...
                    param)
                found = 1;
                break
            end
        end
        if found
            continue
        end
    end

    newchild = [];
   if ~isequal(param,child.parameter,options)
        [newchild, found] = child.generalize(param,root,0,options);
        if isempty(newchild)
            continue
        end

        newchild.memory = child.memory;
        if isempty(found)
            newchild.occurrences = child.occurrences; % This might be avoided in order to get specific classes
        else
            found
        end
        child = newchild;
    end

    for j = 1:length(to.occurrences)
        if isequal(to.occurrences(j).pattern,child)
            stop = 1;
            break
        end
    end
    if stop
        continue
    end

    occ2 = child.occurrence(occ,succ);
    if pat.verbose
        display(occ2);
    end

    if ~isempty(newchild)
        abstract = 0;
        succto = succ;
        if isa(succ,'pat.syntagm')
            succto = succ.to;
        end
        for j = 1:length(succto.occurrences)
            if isequal(succto.occurrences(j),occ2)
                continue
            end
            if ~more_specific_suffix(occ2,succto.occurrences(j))
                child.abstract = true;
                abstract = 1;
                break
            end
        end
        if ~abstract
            child
        end
    end
end


%%
function patterns = learn(memoparam,param,paramstruct,occ,succ,objpat,root,search,patterns)
[idx, memo, value] = memoparam.find(param,[]);
if ~isempty(value)
    if search && ~isempty(memo)
        if iscell(memo)
            if (isa(memo{1}{2},'pat.syntagm') && ...
                    ~isequal(succ.to,memo{1}{2}.to) && ...
                    ~isequal(memo{1}{1},occ)) || ...
                    (isa(memo{1}{2},'pat.event') && ...
                    ~isequal(succ,memo{1}{2}))
                % If we find a new parametric repetition,
                % we can check that the pattern is valid (not redundant)

                param = succ.parameter;
                options = [];
                for i = 1:length(memo)
                    param = memo{i}{2}.parameter.common(param,options);
                end
                param.fields{1}.value = []; %%%%%%%%%%%%%%%%%%%%%%%%%

                if isempty(param.fields{2}.value) || ...
                        (~isequal(objpat,root) && ...
                         isempty(param.fields{1}.inter.value))
                    return
                end

                % Pattern closure check
                nboccs = length(memo) + 1;
                if closuretest(occ,succ,param,objpat,nboccs)
                    % Success! New pattern created as an extension of the parent
 
                    for i = 1:length(objpat.children)
                        childiparam = objpat.children{i}.parameter;
                        if isequal(childiparam.fields{1}.inter,param.fields{1}.inter) && ...
                                isequal(childiparam.fields{2}.value,param.fields{2}.value)
                            return
                        end
                    end

                    child = pat.pattern(root,objpat,param,objpat.memory);

                    patterns(end+1) = child;

                    % Memorising all possible continuations of extended pattern
                    % occurrences
                    for i = 1:length(memo)
                        old = child.occurrence(memo{i}{1},memo{i}{2},1);

                        last1 = memo{i}{2};
                        if ~isa(last1,'pat.event')
                            last1 = last1.to;
                        end
                        for j = 1:length(last1.from)
                            oldsucc = last1.from(j);
                            patterns = learn(child.memory.fields{1}.inter,...
                                oldsucc.parameter.fields{1}.inter,paramstruct,...
                                old,oldsucc,child,root,false,patterns);
                            patterns = learn(child.memory.fields{2},...
                                oldsucc.parameter.fields{2},paramstruct,...
                                old,oldsucc,child,root,false,patterns);
                        end
                    end

                    % Creating the new pattern occurrence
                    occ2 = child.occurrence(occ,succ,1);

                    succto = succ;
                    if isa(succ,'pat.syntagm')
                        succto = succ.to;
                    end
                    abstract = 0;
                    for j = 1:length(succto.occurrences)
                        if isequal(succto.occurrences(j),occ2)
                            continue
                        end
                        if ~more_specific_suffix(occ2,succto.occurrences(j))
                            child.abstract = true;
                            abstract = 1;
                            break
                        end
                    end

                    if ~abstract
                        child
                    end
                end

                % Storing the new continuation in the memory table
                if isempty(idx)
                    memoparam.values(end+1) = value;
                    memoparam.content{end+1} = {{occ,succ}};
                elseif iscell(memoparam.content{idx}) && ...
                        ~isequal(succ,...
                        memoparam.content{idx}{end}{2})
                    memoparam.content{idx}{end+1} = {occ,succ};
                end
            end
        elseif isa(memo,'seq.paramstruct') ...
                && ~succ.parameter.implies(memo)
            found = 0;
            for i = 1:length(parent.children)
                if succ.parameter.implies(...
                        parent.children{i}.parameter)
                    found = 1;
                    break
                end
            end
            if ~found
                newparam = memo.common(succ.parameter);
                if newparam.isdefined(parent)
                    '>>>>>'
                    newpat = pat.pattern(root,parent,newparam,...
                        parent.memory)

                    patterns(end+1) = newpat;

                    if isempty(idx)
                        memoparam.values(end+1) = value;
                        memoparam.content{end+1} = newpat.parameter;
                    else
                        memoparam.content{idx} = newpat.parameter;
                    end
                    newpat.occurrence(occ,succ)
                end
            end
        end
    elseif isempty(idx)
        memoparam.values(end+1) = value;
        memoparam.content{end+1} = {{occ,succ}};
    elseif iscell(memoparam.content{idx}) && ...
            ~isequal(succ,...
            memoparam.content{idx}{end}{2})
        memoparam.content{idx}{end+1} = {occ,succ};
    end
end

%%
function test = closuretest(pref,suff,param,parent,nboccs)
test = 1;

note = suff;
if isa(note,'pat.syntagm')
    note = note.to;
end
if pref.pattern.implies(parent)
    for i = 1:length(note.occurrences)
        if note.occurrences(i).parameter.implies(param)
            prefi = note.occurrences(i).prefix;
            if prefi.pattern.implies(parent,prefi,pref) && ...
                    length(note.occurrences(i).pattern.occurrences) >= nboccs
                test = 0;
                return
            end
        end
    end
end

%%
function res = more_specific_suffix(occ2,occ1)
res = (~isempty(occ2.parameter.fields{1}.inter) && ...
    ~isempty(occ2.parameter.fields{1}.inter.value) && ...
    (isempty(occ1.parameter.fields{1}.inter) || ...
    isempty(occ1.parameter.fields{1}.inter.value))) || ...
    (~isempty(occ2.parameter.fields{2}.value) && ...
    isempty(occ1.parameter.fields{2}.value));

%%
function res = more_specific(occ1,occ2)
if isempty(occ2.prefix)
    res = true;
    return
elseif isempty(occ1.prefix)
    res = false;
    return
end
if (~isempty(occ2.parameter.fields{1}.inter) && ...
        ~isempty(occ2.parameter.fields{1}.inter.value) && ...
        (isempty(occ1.parameter.fields{1}.inter) || ...
        isempty(occ1.parameter.fields{1}.inter.value))) || ...
        (~isempty(occ2.parameter.fields{2}.value) && ...
        isempty(occ1.parameter.fields{2}.value))
    res = false;
    return
end
res = more_specific(occ1.prefix,occ2.prefix);

%%
function [p,s] = defaultParse(s,ps)
% The parametrical description of the new event
p = ps.type2val; % Defined from the parametrical space ps
% ... by associating a value to each parameter (here just one, the
% "dimension")
p = p.setfield('dimension',seq.paramval(ps.getfield('dimension'),s));