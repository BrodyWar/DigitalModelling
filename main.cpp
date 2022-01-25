#define _USE_MATH_DEFINES

#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include <cmath>

double g = 9.81, k = 1, r = 0.6, len = 15;

bool f1 = false, f2 = false;

struct State {
    double x, y;
    State() : x(0), y(0) {};
    State(double x, double y) :x(x), y(y) {}
    friend State operator+(const State& a, const State& b) {
        return State(a.x + b.x, a.y + b.y);
    }
    friend State operator*(const State& a, const State& b) {
        return State(a.x * b.x, a.y * b.y);
    }
    friend State operator-(const State& a, const State& b) {
        return State(a.x - b.x, a.y - b.y);
    }
    friend State operator*(const State& a, const double& b) {
        return State(a.x * b, a.y * b);
    }
    friend State operator*(const double& b, const State& a) {
        return State(a.x * b, a.y * b);
    }
};


struct Shape {
    double x[4], y[4];
    double a, b;
    Shape(int a, int b) :a(a), b(b) {
        x[0] = x[1] = -a / 2;
        x[2] = x[3] = a / 2;
        y[0] = y[3] = -b / 2;
        y[1] = y[2] = b / 2;
    }
    Shape() {
        a = 40, b = 20;
        x[0] = x[1] = 0;
        x[2] = x[3] = 0;
        y[0] = y[3] = 0;
        y[1] = y[2] = 0;
    }

    friend Shape operator+(const Shape& a, const Shape& b) {
        Shape c;
        for (int i = 0; i < 4; i++) {
            c.x[i] = a.x[i] + b.x[i];
            c.y[i] = a.y[i] + b.y[i];
        }
        return c;
    }

    friend Shape operator-(const Shape& a, const Shape& b) {
        Shape c;
        for (int i = 0; i < 4; i++) {
            c.x[i] = a.x[i] - b.x[i];
            c.y[i] = a.y[i] - b.y[i];
        }
        return c;
    }

    friend Shape operator*(const Shape& a, const Shape& b) {
        Shape c;
        for (int i = 0; i < 4; i++) {
            c.x[i] = a.x[i] * b.x[i];
            c.y[i] = a.y[i] * b.y[i];
        }
        return c;
    }

    friend Shape operator*(const Shape& a, const double& b) {
        Shape c;
        for (int i = 0; i < 4; i++) {
            c.x[i] = a.x[i] * b;
            c.y[i] = a.y[i] * b;
        }
        return c;
    }
    friend Shape operator*(const double& a, const Shape& b) {
        Shape c;
        for (int i = 0; i < 4; i++) {
            c.x[i] = a * b.x[i];
            c.y[i] = a * b.y[i];
        }
        return c;
    }
};

struct Wheels {
    State pos[2],mov[2];
    Wheels() {
        pos[0] = pos[1] = mov[0] = mov[1] = State();
    }

    Wheels(double a, double b) {
        pos[0] = State(-a,-b);
        pos[1] = State(a, -b);
    }

    friend Wheels operator+(const Wheels& a, const Wheels& b) {
        Wheels c;
        for (int i = 0; i < 2; i++) {
            c.pos[i] = a.pos[i] + b.pos[i];
            c.mov[i] = a.mov[i] + b.mov[i];
        }
        return c;
    }

    friend Wheels operator-(const Wheels& a, const Wheels& b) {
        Wheels c;
        for (int i = 0; i < 2; i++) {
            c.pos[i] = a.pos[i] - b.pos[i];
            c.mov[i] = a.mov[i] - b.mov[i];
        }
        return c;
    }

    friend Wheels operator*(const Wheels& a, const Wheels& b) {
        Wheels c;
        for (int i = 0; i < 2; i++) {
            c.pos[i] = a.pos[i] * b.pos[i];
            c.mov[i] = a.mov[i] * b.mov[i];
        }
        return c;
    }

    friend Wheels operator*(const Wheels& a, const double& b) {
        Wheels c;
        for (int i = 0; i < 2; i++) {
            c.pos[i] = a.pos[i] * b;
            c.mov[i] = a.mov[i] * b;
        }
        return c;
    }

    friend Wheels operator*(const double& b, const Wheels& a) {
        Wheels c;
        for (int i = 0; i < 2; i++) {
            c.pos[i] = a.pos[i] * b;
            c.mov[i] = a.mov[i] * b;
        }
        return c;
    }
};



class Body {
public:
    double m;
    State pos, mov;
    Body(State pos, State v) : pos(pos), mov(v) {
        init();
    }
    Body(State pos) : pos(pos), mov(State(0, 0)) {
        init();
    }
    Body() :Body(State(0, 0), State(0, 0)) {
        init();
    }
    void init() {
        m = 1;
    }
    friend Body operator+(const Body& a, const Body& b) {
        return Body(a.pos + b.pos, b.mov + a.mov);
    }
    friend Body operator*(const Body& a, const Body& b) {
        return Body(a.pos * b.pos, b.mov * a.mov);
    }
    friend Body operator-(const Body& a, const Body& b) {
        return Body(a.pos - b.pos, b.mov - a.mov);
    }
    friend Body operator*(const Body& a, const double& b) {
        return Body(a.pos * b, a.mov * b);
    }
    friend Body operator*(const double& b, const Body& a) {
        return Body(a.pos * b, a.mov * b);
    }
};

class BodyA : public Body {
public:
    Shape *_shape;
    Wheels *_wheels;
    double R[2][2];
    double a, w;
    BodyA(State pos, State v, int ang) : Body(pos, v) {
        init(ang);
 
    }
    BodyA(State pos, int ang) : Body(pos) {
        init(ang);

    }
    BodyA() :Body(State(0, 0), State(0, 0)) {
        a = 0;
        w = 0.0;
        _shape = new Shape();
        _wheels = new Wheels();
    }
    void init(int ang) { 
        a = ang;
        w = 0.0;
        _shape = new Shape(40, 20);
        _wheels = new Wheels(20, 30);

    }

    friend BodyA operator+(const BodyA& a, const BodyA& b) {
        BodyA h;
        h.pos = a.pos + b.pos;
        h.mov = b.mov + a.mov;
        h.a = a.a + b.a;
        h.w = a.w + b.w;
        *h._shape = *a._shape + *b._shape;
        *h._wheels = *a._wheels + *b._wheels;
        return h;
    }
    friend BodyA operator*(const BodyA& a, const BodyA& b) {
        BodyA h;
        h.pos = a.pos * b.pos;
        h.mov = b.mov * a.mov;
        h.a = a.a * b.a;
        h.w = a.w * b.w;
        *h._shape = *a._shape * *b._shape;
        *h._wheels = *a._wheels * *b._wheels;
        return h;
    }
    friend BodyA operator-(const BodyA& a, const BodyA& b) {
        BodyA h;
        h.pos = a.pos - b.pos;
        h.mov = b.mov - a.mov;
        h.w = a.w - b.w;
        h.a = a.a - b.a;
        *h._shape = *a._shape - *b._shape;
        *h._wheels = *a._wheels - *b._wheels;
        return h;
    }
    friend BodyA operator*(const BodyA& a, const double& b) {
        BodyA h;
        h.pos = a.pos * b;
        h.mov = b * a.mov;
        h.w = a.w * b;
        h.a = a.a * b;
        *h._shape = *a._shape * b;
        *h._wheels = *a._wheels * b;
        return h;
    }
    friend BodyA operator*(const double& b, const BodyA& a) {
        BodyA h;
        h.pos = a.pos * b;
        h.mov = b * a.mov;
        h.w = a.w * b;
        h.a = a.a * b;
        *h._shape = *a._shape * b;
        *h._wheels = *a._wheels * b;
        return h;
    }
    Shape global(Shape* shape) {
        Shape elem;
        double a = this->a / 180.0 * M_PI;
        for (int i = 0; i < 4; ++i) {
            elem.x[i] = this->pos.x + shape->x[i] * std::cos(a) - shape->y[i] * std::sin(a);
            elem.y[i] = this->pos.y + shape->x[i] * std::sin(a) + shape->y[i] * std::cos(a);
        }
        return elem;
    }

    Wheels global(Wheels* wheels) {
        Wheels body;
        double a = this->a / 180.0 * M_PI;
        for (int i = 0; i < 2; i++) {
            body.pos[i].x = this->pos.x + wheels->pos[i].x * std::cos(a) - wheels->pos[i].y * std::sin(a);
            body.pos[i].y = this->pos.y + wheels->pos[i].x * std::sin(a) + wheels->pos[i].y * std::cos(a);
        }
        return body;
    }
};



void f(BodyA* body, BodyA* dot) {
    double m = body[0].m;
    double moment = (1.0 / 12.0) * (std::pow(body[0]._shape->a, 2) + std::pow(body[0]._shape->b, 2)) * m;

    Shape el = body[0].global(body[0]._shape);
    Wheels wheels = body[0].global(body[0]._wheels);
    double a = body[0].a / 180.0 * M_PI; // translate grad to rad


    


    dot[0].a = body[0].w * 180 / M_PI; // translate rad to grad

    dot[0].pos.x = body[0].mov.x;
    dot[0].pos.y = body[0].mov.y;
    dot[0].mov.x = 0;

    dot[0].mov.y =- g -k / m * (body[0]._shape->y[0] - body[0]._wheels->pos[0].y - len) - double(r * body[0].mov.y)
        - k / m * (body[0]._shape->y[3] - body[0]._wheels->pos[1].y - len) -double(r * body[0].mov.y);

    dot[0]._wheels->pos[0].x = 0;
    dot[0]._wheels->pos[0].y = body[0]._wheels->mov[0].y;
    dot[0]._wheels->mov[0].x = 0;
    dot[0]._wheels->mov[0].y = - g + k / m * (body[0]._shape->y[0] - body[0]._wheels->pos[0].y - len) +double(r * body[0].mov.y);


    dot[0]._wheels->pos[1].x = 0;
    dot[0]._wheels->pos[1].y = body[0]._wheels->mov[1].y;
    dot[0]._wheels->mov[1].x = 0;
    dot[0]._wheels->mov[1].y = - g + k / m * (body[0]._shape->y[3] - body[0]._wheels->pos[1].y - len) +double(r * body[0].mov.y);


}

void step(BodyA st[], BodyA* body, BodyA* k, double h) {
    for (int i = 0; i < 1; ++i)
        st[i] = body[i] + h / 2 * (k[i]);
}

void RungeKutta(BodyA* body, double h, void (*f)(BodyA*, BodyA*)) {
    BodyA k1[1], k2[1], k3[1], k4[1], st[1];
    f(body, k1);
    step(st, body, k1, h);
    f(st, k2);

    step(st, body, k2, h);
    f(st, k3);

    step(st, body, k3, h * 2);
    f(st, k4);

    for (int i = 0; i < 1; ++i)
        body[i] = body[i] + (h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
}


int main()
{
    GLFWwindow* window;


    if (!glfwInit())
        return -1;

    window = glfwCreateWindow(800, 640, "Hello World", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    BodyA body[1];
    body[0] = BodyA(State(0, 50), 40); // center (x,y), angle[grad]
    double h = 0.015;

    glfwMakeContextCurrent(window);
    glScaled(0.01f, 0.01f, 0.0);
    while (!glfwWindowShouldClose(window))
    {
        
        int state = glfwGetKey(window, GLFW_KEY_SPACE);

        if (state == GLFW_PRESS)
        {
            body[0] = BodyA(State(0, 50), 0);
        }

        if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
        {
            body[0].w += 0.3;
        }

        else if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
        {
            body[0].w -= 0.3;
        }
        else
            body[0].w += 0;

        if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
        {
            body[0].mov.x = 10;
        }

        else if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
        {
            body[0].mov.x = -10;
        }
        else
            body[0].mov.x = 0;


        {
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glClearColor(0.2f, 0.3f, 0.6f, 0.0f);


            Shape el = body[0].global(body[0]._shape);
            Wheels wheels = body[0].global(body[0]._wheels);

            glLineWidth(3);
            glBegin(GL_QUADS);
            glColor3f(0.0f, 0.0f, 0.0f);

            glVertex2f(el.x[0], el.y[0]);
            glVertex2f(el.x[1], el.y[1]);
            glVertex2f(el.x[2], el.y[2]);
            glVertex2f(el.x[3], el.y[3]);
            glEnd();

            glLineWidth(1);
            glBegin(GL_LINES);
            glColor3f(1.0f, 1.0f, 1.0f);
            glVertex2f(-70, 0);
            glVertex2f(70, 0);
            glColor3f(0.4f, 1.0f, 0.7f);
            glVertex2f(wheels.pos[0].x, wheels.pos[0].y);
            glVertex2f(el.x[0], el.y[0]);

            glVertex2f(wheels.pos[1].x, wheels.pos[1].y);
            glVertex2f(el.x[3], el.y[3]);

            glEnd();


            glPointSize(4);
            glBegin(GL_POINTS);
            glColor3f(0.0f, 0.0f, 0.0f);
            glVertex2f(wheels.pos[0].x, wheels.pos[0].y);
            glVertex2f(wheels.pos[1].x, wheels.pos[1].y);
            glEnd();


            glfwSwapBuffers(window);
            glFlush();
            glfwPollEvents();
        }



        double m = body[0].m;
        double a = body[0].a / 180.0 * M_PI; // translate grad to rad
        double moment = (1.0 / 12.0) * (std::pow(body[0]._shape->a, 2) + std::pow(body[0]._shape->b, 2)) * m;
        Shape el1 = body[0].global(body[0]._shape);
        Wheels wheels1 = body[0].global(body[0]._wheels);

        double tmp = body[0].w;
        if ((el1.y[0] <= 0 || el1.y[1] <= 0 || el1.y[2] <= 0 || el1.y[3] <= 0) && body[0].mov.y < 0) {
            body[0].mov.y = 0;  
            if (el1.y[2] <= 0)
                body[0].w += -g * (el1.x[2] - body[0].pos.x) / (moment);
            if (el1.y[3] <= 0)
                body[0].w += g * (el1.x[3] - body[0].pos.x) / (moment);
            if (el1.y[1] <= 0)
                body[0].w += -g * (el1.x[1] - body[0].pos.x) / (moment);
            if (el1.y[0] <= 0)
                body[0].w += g * (el1.x[0] - body[0].pos.x) / (moment);
            if (std::fabs(body[0].w - tmp) < 0.1) {
                body[0].w = 0;
            }
        }
        
        double curlen1 = body[0]._shape->y[0] - body[0]._wheels->pos[0].y, curlen2 = body[0]._shape->y[3] - body[0]._wheels->pos[1].y;

        if (wheels1.pos[0].y <= 0) { // left wheel
            if (std::fabs(body[0]._wheels->mov[0].y) < 0.25)
                body[0]._wheels->mov[0].y = 0;

            body[0]._wheels->mov[0].y = -body[0].mov.y -wheels1.pos[0].y;
            body[0].w += k * std::fabs(std::sqrt(std::pow(wheels1.pos[0].y - el1.y[0], 2) + std::pow(wheels1.pos[0].x - el1.x[0], 2)) - len )  * (el1.x[0] - body[0].pos.x) / (moment);

 


        }

        if (wheels1.pos[1].y <= 0) { // right wheel
            if (std::fabs(body[0]._wheels->mov[1].y) < 0.25)
                body[0]._wheels->mov[1].y = 0;

            body[0]._wheels->mov[1].y = -body[0].mov.y -wheels1.pos[1].y;
            body[0].w += k * std::fabs(std::sqrt(std::pow(wheels1.pos[1].y - el1.y[3], 2) + std::pow(wheels1.pos[1].x - el1.x[3], 2)) - len)  * (el1.x[3] - body[0].pos.x) / (moment);



        }
        body[0].w = body[0].w * 0.8;

        RungeKutta(body, h, f);
    }
    glfwTerminate();
    return 0;
}
