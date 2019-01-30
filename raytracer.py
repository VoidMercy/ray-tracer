import copy, math

# TODO: add index of refraction
# add fresnel's equation
# add physical sun as a white sphere

bg_color = [153, 194, 255]

class Vector:
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z

	def distance(self, point2):
		return math.sqrt((self.x - point2.x) ** 2 + (self.y - point2.y) ** 2 + (self.z - point2.z) ** 2)

	def __sub__(self, point2):
		return Vector(self.x - point2.x, self.y - point2.y, self.z - point2.z)

	def dot(self, point2):
		return self.x * point2.x + self.y * point2.y + self.z * point2.z

	def dotvec(self, vec2):
		return Vector(self.x * vec2.x, self.y * vec2.y, self.z * vec2.z)

	def mag(self):
		return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

	def cross(self, vec2):
		return Vector(self.y * vec2.z - self.z * vec2.y, self.z * vec2.x - self.x * vec2.z, self.x * vec2.y - self.y * vec2.x)

	def __str__(self):
		return "(%s, %s, %s)" % (self.x, self.y, self.z)

	def __div__(self, val):
		return Vector(self.x / val, self.y / val, self.z / val)

	def __add__(self, vec2):
		return Vector(self.x + vec2.x, self.y + vec2.y, self.z + vec2.z)

	def __mul__(self, val):
		return Vector(self.x * val, self.y * val, self.z * val)

	def norm(self):
		return self / self.mag()

	def toarray(self):
		return [self.x, self.y, self.z]

class Line:
	def __init__(self, point1, point2):
		self.point1 = copy.deepcopy(point1)
		self.point2 = copy.deepcopy(point2)
		self.direction = Vector(point2.x - point1.x, point2.y - point1.y, point2.z - point1.z)

class Sphere:

	def __init__(self, location, r, reflectivity, transparency, color):
		self.location = location
		self.r = r
		self.reflectivity = reflectivity
		self.transparency = transparency
		self.color = color

	def intersect(self, line):
		A = line.direction.x ** 2 + line.direction.y ** 2 + line.direction.z ** 2
		B = 2 * line.direction.x * (line.point1.x - self.location.x) + 2 * line.direction.y * (line.point1.y - self.location.y) + 2 * line.direction.z * (line.point1.z - self.location.z)
		C = line.point1.x ** 2 + line.point1.y ** 2 + line.point1.z ** 2 - 2 * line.point1.x * self.location.x - 2 * line.point1.y * self.location.y - 2 * line.point1.z * self.location.z - self.r**2 + self.location.x ** 2 + self.location.y ** 2 + self.location.z ** 2

		temp = B**2 - 4 * A * C
		if temp < 0:
			return []
		t1 = (-B + math.sqrt(temp)) / (2. * A)
		t2 = (-B - math.sqrt(temp)) / (2. * A)
		ret = []
		if t1 > 0.00000000001:
			p = Vector(line.direction.x * t1 + line.point1.x, line.direction.y * t1 + line.point1.y, line.direction.z * t1 + line.point1.z)
			ret.append(p)
		if t2 > 0.00000000001:
			p = Vector(line.direction.x * t2 + line.point1.x, line.direction.y * t2 + line.point1.y, line.direction.z * t2 + line.point1.z)
			ret.append(p)

		return ret

	def get_normal(self, point):
		return Vector(2 * (point.x - self.location.x), 2 * (point.y - self.location.y), 2 * (point.z - self.location.z))

	def get_color(self, closest_point, lightsource, objects, dirvec, depth):
		pixel_value = Vector(0, 0, 0)
		if (self.reflectivity == 0 and self.transparency == 0) or depth > 5:
			path = Line(closest_point, lightsource)
			clear = True
			for obj in objects:
				if obj != self:
					temp_points = obj.intersect(path)
					if len(temp_points) > 0:
						clear = False
						break

			if clear:

				vec1 = self.get_normal(closest_point)
				vec2 = lightsource - closest_point
				multiplier = max(vec1.dot(vec2) / (vec1.mag() * vec2.mag()), 0)
			else:
				multiplier = 0

			toadd = Vector(closest_object.color[0], closest_object.color[1], closest_object.color[2])

			percentage = 0.75

			pixel_value += toadd * multiplier * percentage + toadd * (1 - percentage)
		else:

			reflection = Vector(0, 0, 0)
			refraction = Vector(0, 0, 0)
			inside = False
			if dirvec.dot(self.get_normal(closest_point)) > 0:
				inside = True
			if self.reflectivity != 0:
				N = self.get_normal(closest_point).norm()
				if inside:
					N *= -1
				L = dirvec
				R = L - N * (2 * N.dot(L))
				ray = Line(closest_point, closest_point + R)

				second_intersection = None
				closest_distance = 99999999999999
				second_object = None

				# get the closest intersection point between the ray and the scene
				for obj in objects:
					temp_points = obj.intersect(ray)
					for i in temp_points:
						dist = closest_point.distance(i)
						if dist < closest_distance:
							closest_distance = dist
							second_intersection = i
							second_object = obj
				if second_intersection != None:
					reflection = second_object.get_color(second_intersection, lightsource, objects, second_intersection - closest_point, depth + 1) * self.reflectivity
				else:
					reflection = Vector(bg_color[0], bg_color[1], bg_color[2]) * self.reflectivity
			if self.transparency != 0:
				ray = Line(closest_point, dirvec + closest_point)
				second_intersection = None
				closest_distance = 99999999999999
				second_object = None

				# get the closest intersection point between the ray and the scene
				for obj in objects:
					temp_points = obj.intersect(ray)
					for i in temp_points:
						dist = closest_point.distance(i)
						if dist < closest_distance:
							closest_distance = dist
							second_intersection = i
							second_object = obj
				if second_intersection != None:
					refraction = second_object.get_color(second_intersection, lightsource, objects, second_intersection - closest_point, depth + 1) * self.transparency
				else:
					refraction = Vector(bg_color[0], bg_color[1], bg_color[2]) * self.transparency

			pixel_value = (reflection * 0.5 + refraction * 0.5).dotvec(Vector(self.color[0], self.color[1], self.color[2]) / 255.)

		return pixel_value


class Plane:

	def __init__(self, normal, point, color):
		self.normal = normal
		self.point = point
		self.color = color

	def intersect(self, line):
		temp = (self.normal.x * line.direction.x + self.normal.y * line.direction.y + self.normal.z * line.direction.z)
		if temp == 0:
			return []
		t = (self.normal.x * (self.point.x - line.point1.x) + self.normal.y * (self.point.y - line.point1.y) + self.normal.z * (self.point.z - line.point1.z)) / float(temp)
		if t <= 0:
			return []
		return [Vector(line.direction.x * t + line.point1.x, line.direction.y * t + line.point1.y, line.direction.z * t + line.point1.z)]

	def get_normal(self, point):
		return self.normal

	def get_color(self, closest_point, lightsource, objects, dirvec, depth):
		color = (int(round(int(closest_point.z) / 1000.)) + int(round(int(closest_point.x) / 1000.))) % 2
		if color & 1 == 1:
			color = Vector(255, 70, 70)
		else:
			color = Vector(0, 230, 0)
		path = Line(closest_point, lightsource)
		clear = True
		for obj in objects:
			if obj != self:
				temp_points = obj.intersect(path)
				if len(temp_points) > 0:
					clear = False
					break

		if clear:
			vec1 = self.get_normal(closest_point)
			vec2 = lightsource - closest_point
			multiplier = max(vec1.dot(vec2) / (vec1.mag() * vec2.mag()), 0)
		else:
			multiplier = 0

		toadd = color

		percentage = 0.70

		pixel_value = toadd * multiplier * percentage + toadd * (1 - percentage)
		return pixel_value

objects = []
ground = Plane(Vector(0., 1., 0.), Vector(0, 0, 0), [255, 0, 0])
sphere = Sphere(Vector(0., 500., -300.), 100., 1, 0, [255, 255, 255])
sphere2 = Sphere(Vector(200., 500., -300.), 100., 1, 0, [255, 255, 255])
sphere3 = Sphere(Vector(100., 570., -300.), 40., 1, 0, [255, 255, 255])
objects.append(ground)
objects.append(sphere)
objects.append(sphere2)
# objects.append(sphere3)

# point of view
viewpoint = Vector(-300, 700, 700.)
projection_point = Vector(-180, 560, 100.)
d_projection = viewpoint - projection_point

inc_x = ground.normal.cross(d_projection).norm()
inc_y = d_projection.cross(inc_x).norm()

# source of light
lightsource = Vector(2500., 2000., 3000.)

width = 1024
height = 768
step_size_x = 1
step_size_y = 1

out = open("image.ppm", "w")

out.write("P6 %d %d %d " % (width, height, 255))

for y in range(height):
	for x in range(width):

		pixel_value = [0, 0, 0]

		# create a ray from the viewpoint to the pixel
		cur_projection = (inc_x * step_size_x * (-(width/2.0) + x)) + (inc_y * step_size_y * (height - (height/2.0) - y)) + projection_point
		ray = Line(viewpoint, cur_projection)
		closest_point = None
		closest_distance = 99999999999999
		closest_object = None

		# get the closest intersection point between the ray and the scene
		for obj in objects:
			temp_points = obj.intersect(ray)
			for i in temp_points:
				dist = viewpoint.distance(i)
				if dist < closest_distance:
					closest_distance = dist
					closest_point = i
					closest_object = obj

		if closest_point != None:

			pixel_value = closest_object.get_color(closest_point, lightsource, objects, closest_point - viewpoint, 0).toarray()
		else:
			pixel_value = bg_color

		out.write(''.join([chr(int(i)) for i in pixel_value]))

out.close()


'''

Center: (x0, y0, z0)
Radius: r
Sphere: (x-x0)^2 + (y-y0)^2 + (z-z0)^2 = r^2

Gradient: <2*(x - x0), 2*(y-y0), 2*(z-z0)>

Line:
P1: (x1, y1, z1)
P2: (x2, y2, z2)
D: (dx, dy, dz) = (x2 - x1, y2 - y1, z2 - z1)

x = dx * t + x1
y = dy * t + y1
z = dz * t + z1

x^2 - (2 * x * x0) + x0^2 + y^2 - (2 * y * y0) + y0^2 + z^2 - (2 * z * z0) + z0^2 = r^2
x^2 - 2*x*x0 + y^2 - 2*y*y0 + z^2 - 2*z*z0 = r^2 - x0^2 - y0^2 - z0^2

dx^2 * t^2 + 2*x1*dx*t + x1^2 - 2*x0*dx*t - 2*x0*x1
t^2 * (dx^2 + dy^2 + dz^2) + t * (2*x1*dx - 2*x0*dx + 2*y1*dy - 2*y0*dy + 2*z1*dz - 2*z0*dz) + x1^2 + y1^2 + z1^2 - 2*x0*x1 - 2*y0*y1 - 2*z0*z1 = r^2 + x0^2 + y0^2 + z0^2

A = dx^2 + dy^2 + dz^2
B = 2*x1*dx - 2*x0*dx + 2*y1*dy - 2*y0*dy + 2*z1*dz - 2*z0*dz
  = 2 * dx * (x1 - x0) + 2 * dy * (y1 - y0) + 2 * dz * (z1 - z0)
C = x1^2 + y1^2 + z1^2 - 2*x0*x1 - 2*y0*y1 - 2*z0*z1 - r^2 + x0^2 + y0^2 + z0^2


Plane:
Normal: (nx, ny, nz)
Point: (x0, y0, z0)
nx(x - x0) + ny(y - y0) + nz(z - z0) = 0

Line:
P1: (x1, y1, z1)
P2: (x2, y2, z2)
D: (dx, dy, dz) = (x2 - x1, y2 - y1, z2 - z1)

x = dx * t + x1
y = dy * t + y1
z = dz * t + z1

nx * x - nx * x0 + ny * y - ny * y0 + nz * z - nz * z0 = 0
nx * (dx * t + x1) + ny * (dy * t + y1) + nz * (dz * t + z1) - nx * x0 - ny * y0 - nz * z0 = 0

t * (nx * dx + ny * dy + nz * dz) + nx * x1 + ny * y1 + nz * z1 - nx * x0 - ny * y0 - nz * z0 = 0

t = (nx * (x0 - x1) + ny * (y0 - y1) + nz * (z0 - z1)) / (nx * dx + ny * dy + nz * dz)

Viewpoint: (vx, vy, vz)
Projection center: (px, py, pz)

Direction vector: (px - vx, py - vy, pz - vz) = (dx, dy, dz)
Projection plane: dx(x - px) + dy(y - py) + dz(z - pz) = 0

Ground normal vector: <0, 1, 0>

Vector1:
x = x1 * t + px
y = y1 * t + py
z = z1 * t + pz

(x1, y1, z1) = direction vector X ground normal vector
'''