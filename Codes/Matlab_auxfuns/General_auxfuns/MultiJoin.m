function Joined = MultiJoin(varargin,keys)
    Joined = varargin{1};
    n = length(varargin);
    for k = 2:n
        Joined = join(Joined, varargin{k},'Keys',keys);
    end
end