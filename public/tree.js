let treeData = {
    name: "But initial",
    children: []
  };
  
  const width = 800;
  const height = 600;
  
  const svg = d3.select("#tree")
    .append("svg")
    .attr("width", width)
    .attr("height", height)
    .append("g")
    .attr("transform", "translate(50, 50)");
  
  function update(treeData) {
    svg.selectAll("*").remove();
  
    const treeLayout = d3.tree().size([width - 100, height - 100]);
    const root = d3.hierarchy(treeData);
    treeLayout(root);
  
    svg.selectAll('.link')
      .data(root.links())
      .enter()
      .append('line')
      .classed('link', true)
      .attr('x1', d => d.source.x)
      .attr('y1', d => d.source.y)
      .attr('x2', d => d.target.x)
      .attr('y2', d => d.target.y)
      .style('stroke', '#555');
  
    const node = svg.selectAll('.node')
      .data(root.descendants())
      .enter()
      .append('g')
      .classed('node', true)
      .attr('transform', d => `translate(${d.x},${d.y})`);
  
    node.append('circle')
      .attr('r', 5)
      .style('fill', d => d.data.success === 'full' ? 'green' : 'orange');
  
    // Ajouter un carré pour les succès partiels
    node.filter(d => d.data.success === 'partial')
      .append('rect')
      .attr('width', 10)
      .attr('height', 10)
      .attr('x', -5)
      .attr('y', -5)
      .style('fill', 'orange');
  
    node.append('text')
      .attr('dy', -10)
      .attr('text-anchor', 'middle')
      .text(d => `${d.data.name} (${d.data.rule || 'But initial'})`);
  }
  
  function updateTreeWithRules(but, rule, status) {
    // Mettre à jour l'arbre avec les résultats de la trace
    // Créer un nouvel objet représentant le nouveau nœud de l'arbre
    const newNode = {
      name: but,
      rule: rule,
      success: status,
      children: []
    };
  
    // Trouver le nœud parent où attacher le nouveau nœud
    const parent = findParentNode(treeData, rule);
  
    // Attacher le nouveau nœud au parent
    parent.children.push(newNode);
  
    // Mettre à jour l'arbre avec le nouveau nœud
    update(treeData);
  }
  
  function findParentNode(node, rule) {
    // Parcourir l'arbre de manière récursive pour trouver le parent du nœud actuel
    if (node.rule === rule) {
      return node;
    } else {
      for (const child of node.children) {
        const parent = findParentNode(child, rule);
        if (parent) {
          return parent;
        }
      }
    }
  }
  
  // Fonction pour mettre à jour l'arbre initial avec le but initial
  update(treeData);
  