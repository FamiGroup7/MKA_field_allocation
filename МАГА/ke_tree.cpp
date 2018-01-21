#include "ke_tree.h"

ke_tree::ke_tree()
{
}

ke_tree::~ke_tree()
{
	points = NULL;
	kes = NULL;
}

void ke_tree::init(double x1, double x2, double y1, double y2, double z1, double z2, vector<Point>*points, vector<nvtr>* kes) {
	this->delim = 0;
	this->points = points;
	this->kes = kes;
	this->info = field(x1, x2, y1, y2, z1, z2, -1, -1);
	x1_outer = x1; x2_outer = x2;
	y1_outer = y1; y2_outer = y2;
	z1_outer = z1; z2_outer = z2;
}

void ke_tree::init(int parentDelim, ke_tree &parent, bool first) {
	this->info.x1 = parent.info.x1;
	this->info.x2 = parent.info.x2;
	this->info.y1 = parent.info.y1;
	this->info.y2 = parent.info.y2;
	this->info.z1 = parent.info.z1;
	this->info.z2 = parent.info.z2;
	this->points = parent.points;
	this->kes = parent.kes;
	switch (parentDelim)
	{
	case 0: delim = 1;
		first ?
			this->info.x2 = (parent.info.x1 + parent.info.x2) / 2.0 :
			this->info.x1 = (parent.info.x1 + parent.info.x2) / 2.0;
		break;
	case 1: delim = 2;
		first ?
			this->info.y2 = (parent.info.y1 + parent.info.y2) / 2.0 :
			this->info.y1 = (parent.info.y1 + parent.info.y2) / 2.0;
		break;
	case 2: delim = 0;
		first ?
			this->info.z2 = (parent.info.z1 + parent.info.z2) / 2.0 :
			this->info.z1 = (parent.info.z1 + parent.info.z2) / 2.0;
		break;
	default:
		throw new exception("error construct help struct for ke");
	}
	x1_outer = info.x1; x2_outer = info.x2;
	y1_outer = info.y1; y2_outer = info.y2;
	z1_outer = info.z1; z2_outer = info.z2;
}

Point ke_tree::centerOfKe(int iKe) {
	double xCenter, yCenter, zCenter;
	xCenter = points->at(kes->at(iKe).uzel[0]).x + points->at(kes->at(iKe).uzel[7]).x;
	yCenter = points->at(kes->at(iKe).uzel[0]).y + points->at(kes->at(iKe).uzel[7]).y;
	zCenter = points->at(kes->at(iKe).uzel[0]).z + points->at(kes->at(iKe).uzel[7]).z;
	return Point(xCenter / 2, yCenter / 2, zCenter / 2);
}

bool ke_tree::contains(Point point) {
	return point.inside(info.x1, info.x2, info.y1, info.y2, info.z1, info.z2);
}

bool ke_tree::containsInOuter(Point point) {
	return point.inside(x1_outer, x2_outer, y1_outer, y2_outer, z1_outer, z2_outer);
}

void ke_tree::addKe(int iKe) {
	addKe(iKe, centerOfKe(iKe));
}

int ke_tree::findKe(Point point) {
	if (!containsInOuter(point)) return -1;
	if (child1 == NULL && child2 == NULL) {
		for each (int var in elements)
		{
			Point leftDown = points->at(kes->at(var).uzel[0]);
			Point rightUp = points->at(kes->at(var).uzel[7]);
			if (point.inside(leftDown.x, rightUp.x, leftDown.y, rightUp.y, leftDown.z, rightUp.z)) {
				return var;
			}
		}
		return -1;
	}
	int result = child1->findKe(point);
	if (result == -1) {
		result = child2->findKe(point);
	}
	return result;
}

bool ke_tree::addKe(int iKe, Point point) {
	bool added = false;
	if (child1 == NULL && child2 == NULL) {
		if (!contains(point)) {
			throw new exception("dsadf");
		}
		elements.push_back(iKe);
		added = true;
		if (elements.size() >= HELP_STRUCT_KE_MAX_SIZE) {
			child1 = new ke_tree();
			child2 = new ke_tree();
			child1->init(delim, *this, true);
			child2->init(delim, *this, false);
			for each (int var in elements)
			{
				addKe(var);
			}
			elements.clear();
		}
	}
	else {
		if (child1->contains(point)) {
			added = child1->addKe(iKe, point);
		}
		else if (child2->contains(point)) {
			added = child2->addKe(iKe, point);
		}
		else throw new exception("cannot add iKe");
	}
	if (added) {
		expand(points->at(kes->at(iKe).uzel[0]), points->at(kes->at(iKe).uzel[7]));
	}
	return added;
}

void ke_tree::expandOuter(double probablyNewOuter, double&outer, int comparison) {
	if (probablyNewOuter*comparison < outer*comparison) {
		outer = probablyNewOuter;
	}
}

void ke_tree::expand(Point leftDown, Point rightUp) {
	expandOuter(leftDown.x, x1_outer, 1);
	expandOuter(leftDown.y, y1_outer, 1);
	expandOuter(leftDown.z, z1_outer, 1);
	expandOuter(rightUp.x, x2_outer, -1);
	expandOuter(rightUp.y, y2_outer, -1);
	expandOuter(rightUp.z, z2_outer, -1);
}