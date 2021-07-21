import Two from 'two.js';
import Delaunator from 'delaunator';

let grid = {};

//generate square grid.
let width = 100;
let height = 100;

var two;

export function init(){
  two = new Two({
    fullscreen: true,
    autostart: true
  }).appendTo(document.body);

  two.bind('update', update);
}

export function generate(){
  //generate points in grid
  const spacing = 10;
  grid.points = getGrid(width, height, spacing);

  //triangulate points
  grid.delauney = new Delaunator(grid.points.flat());

  const triangles = grid.delauney.triangles;
  const points = grid.points;
  console.log(points);

  grid.triangle_coordinates = [];
  for (let i = 0; i < triangles.length; i += 3) {
    grid.triangle_coordinates.push([
        points[triangles[i]],
        points[triangles[i + 1]],
        points[triangles[i + 2]]
    ]);
  }

  console.log(grid.triangle_coordinates);

  //voronoi points
}

//generate a grid of points evenly spaced between the width and height.
function getGrid(width, height, spacing){
  const radius = spacing / 2 ;
  let points = [];
  for(let y = radius; y<height; y+=spacing){
    for(let x = radius; x<width; x+=spacing){
      points.push([x + Math.random()*radius - (radius/2),y + Math.random()*radius - (radius/2)]);
    }
  }
  return points;
}

function update(){
  //nothing to update yet.
}

export function draw(){
  grid.triangle_coordinates.forEach(element => {
    two.makeLine(...element[0], ...element[1]);
    two.makeLine(...element[1], ...element[2]);
    two.makeLine(...element[2], ...element[0]);
  });
}






