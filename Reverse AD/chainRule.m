function graph_output = chainRule(graph)
%ChainRule is a function which accumulated derivatives in reverse order, 
% from the to node down to (multiple) start nodes.

% Gustav Hedin, January 2019.

% "graph" is the top-node of the calculation tree. Call this top-node
% "current" for simplicity:
current = graph;
% Set derivative to 1: (df/df = 1)
current.derivative = 1;
% Push derivative to parents, and grandparents etc.. :
push_derivative(current);
% Return modified graph:
graph_output = graph;
end

function push_derivative(node)
% This is a recursive function which pushes partial derivatives to parent
% nodes. (No a breadth-first alrgoritm. Rearrange the nodes to obtain a more
% efficient algorithm.)
    if isempty(node.parents)
        return % If no parents, then we're done!
    end
    % Push derivative from current node to the nodes parents, according to
    % corresponding derivativeOp-function:
    node.derivativeOp(node.derivative, node.parents);
    % Since there exist at least one parent, there exist a "left-side" 
    % parent (node.parents(1)). Call the push_derivative function on the 
    % current node's left-parent.
    push_derivative(node.parents(1))
    % If node has "right-parent", then call the push_derivative function
    % on the current nodes right-parent.
    if size(node.parents,2)==2
        push_derivative(node.parents(2))
    end
    % We have now pused all derivatives down to every start nodes. 
end
