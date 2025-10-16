function M = get_model_by_name(name)
%GET_MODEL_BY_NAME Lookup a distribution definition by its display name
models = dists.registry();
target = string(name);
modelNames = string({models.name});
idx = find(strcmpi(modelNames, target), 1);
if isempty(idx)
    error('Distribution "%s" is not registered.', target);
end
M = models(idx);
end
